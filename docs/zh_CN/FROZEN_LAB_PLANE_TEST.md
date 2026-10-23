# 固定实验室启动面：自由漂移对照

## 启动方式选择

HDF5-v4 固定实验室面输入有以下三种启动方式；`initial_self_field.enabled` 和 `mesh.self_field_feedback` 是与启动方式正交的物理/诊断开关，不是另外两种启动方式。没有显式设置时仍使用原有的束头锚定方式，不改变现有输入卡的默认行为。

| 方式 | `frozen_lab_plane_injection` | `hold_transverse_until_activation` | 启动前轨迹与适用情形 |
|---|---:|---:|---|
| 原有束头锚定 | `false` | 省略或 `false` | 重构共同时间束团，再做原有同步；保留旧卡语义。 |
| 冻结三维（ballistic） | `true` | 省略或 `false` | 按完整三维动量直线传播至启动面，不受网格场力；用作保持原始直线轨迹的对照。 |
| 冻结横向（held） | `true` | `true` | 启动前固定 HDF5 输入面的 x/y，仅按实验室纵向速度推进 z；适合高 boost 回推使大量粒子落到有限横向网格外的工况。 |

对这种高 boost、有限横向孔径的固定实验室面束团，**推荐先选 held**，并将 `activation_plane_z` 设为 `input_plane_z`，以保持输入面 x/y。现有 10 µm、50 A、约 1041 万粒子的 1 m 自由漂移测试中，held 显著减少初始出框，三项网格/CPML 扰动后的出口粒子位置差 RMS 均小于 1 nm，盒内总场能量也未持续失控增长。复现实验包仅在本地保留，不随 APRL 源码仓库发布。这是针对该类数值几何的推荐，**不是**一般物理精确性、出口场或辐射收敛的证明。需要原始三维直线轨道作为参照时仍应跑 ballistic 对照。运行 manifest 的 `resolved_beam_reference.startup_mode` 记录解析后的方式；它不是可填写的输入 YAML 键。

`beam.reference.frozen_lab_plane_injection: true` 是实验性 HDF5-v4 输入模式。可选的 `activation_plane_z` 默认为 `input_plane_z`。在冻结模式中不启用 held 时，以完整三维直线运动把输入面事件同步到计算系时间零；启动面前 x/y/z 都按原动量推进，只关闭电磁力，位于物理区内时仍作为 Maxwell 源。粒子跨过实验室固定启动面的时间步在交点拆开，剩余部分恢复正常 Boris 推进。

新增的可选 `beam.reference.hold_transverse_until_activation: true` 必须与冻结模式一起使用。它在同步及冻结传播中把 x/y 保持为 HDF5 输入面值。冻结纵向速度由实验室 `u_z` 单独决定：把完整动量反变换回实验室，令横向分量为零，再把 `(0, 0, u_{z,lab})` 变换到计算系，按 `v'_{z,hold} = c * u'_{z,hold} / sqrt(1 + u'_{z,hold}²)` 推进 z。相同实验室 `u_z`、不同横向动量给出相同冻结纵向速度。初始 Poisson 的平均速度也采用这一保持模式的场源约定（横向为零、纵向用上述速度）。启动面之后恢复保存的完整三维动量和正常 Boris 推进；过面的一个时间步分段沉积，冻结段不会被错误地当作横向电流，上游粒子探测面会记录冻结运动。因此这是减轻高 boost 回推造成的横向出框的**非物理近似**，启动处实际的横向电流与纵向速度会突然改变，可能产生人为辐射。若需在输入面精确保持 HDF5 x/y，启动面应等于输入面；若设在上游，解冻后的传播可以在到达输入面之前改变 x/y。

例如，在输入面启动横向保持模式的配置片段是：

```yaml
beam:
  reference:
    input_plane_z: -0.05
    frozen_lab_plane_injection: true
    activation_plane_z: -0.05  # 可省略，默认等于 input_plane_z
    hold_transverse_until_activation: true
```

运行自动回归：

```bash
cmake --build build -j4
ctest --test-dir build -R '^required\.frozen-plane\.free-drift' --output-on-failure
ctest --test-dir build -V -R '^required\.frozen-plane\.free-drift$' | rg 'ballistic_|peak_to_initial|max_step_increase_to_initial|high_boost|held_'
```

测试自己生成 4096 条 HDF5-v4 输入面记录，不需要外部束团。原有三维冻结组把启动面放在输入面上游 1 µm；新的三维基线和横向保持组把启动面直接放在输入面。所有组均无磁元件，并设置输入面、其后 12 µm 和 20 µm 的粒子探测面；横向保持组另检查输入面上游 1 µm 的冻结轨道。还设置 12 µm 处的 E/B 场探测面、CPML 和每步场能量统计。原有三组主对照为：

| 组别 | 冻结到启动面 | 初始 Poisson 自场 | 后续 Maxwell 源与反作用 |
|---|---:|---:|---:|
| `frozen-initial-field` | 是 | 开 | 开 |
| `frozen-no-initial-field` | 是 | 关 | 开 |
| `unfrozen-initial-field` | 否 | 开 | 开 |

额外的 `frozen-no-feedback-control` 保留初场与 Maxwell 源、只关闭自场反作用，给出直线轨道基线；`frozen-coupled-auto-z` 检查自动纵向网格。新组包括 `ballistic-at-input`（完整三维基线）、`held-initial-field`（横向保持、初场和反作用开启）、`held-no-initial-field`（横向保持、初场关闭）以及 `held-no-feedback-control`（横向保持、初场开启但无自场反作用）。测试同时在 1 和 2 个 MPI 进程下运行。比较量是探测面记录的 x/y 与输入面动量给出的直线预测的 RMS/最大差，场探测面上的最大 |E| 与 c|B|，以及整个自由漂移过程中物理区场能量的初值、峰值、末值和最大单步增量。回归要求无反作用基线的最大轨道差小于 0.0001 µm，带初场的组别场能量峰值小于初值 2 倍；横向保持＋初场组额外要求峰值小于初值 1.05 倍、最大单步增量小于初值的 5%、场探测面 |E| 与 c|B| 均低于 10⁸ V/m。无初场组从零开始，检查其绝对场能量有限。这些是小样例的故障哨兵，**不是**实际 40 µm 束团的物理收敛判据。

本地 4096 粒子样例中，原有三维无自场反作用组在输入面及两个下游探测面的最大直线轨道差均约 10⁻¹⁴ µm，表明反推与三维漂移在本测试精度内回到输入轨道。原有冻结＋初场组在输入面的最大差约 5.3×10⁻⁷ µm（启动面到输入面之间的自场作用），末端约 5.2×10⁻⁵ µm；场能量峰值/初值约 1.0008。原有冻结但无初场组从零增长到 4.47×10⁻¹¹ J，比带初场组的峰值 2.64×10⁻¹¹ J 高约 69%，表明持续沉积本身**不能**消除零场启动顺态。

同一合成样例的新横向保持组在输入面及上游探测器的最大 x/y 差均为零；`held-initial-field` 在下游 20 µm 处的最大轨道差约 5.05×10⁻⁵ µm，场能量初值 2.566×10⁻¹¹ J、峰值 2.570×10⁻¹¹ J（峰值/初值约 1.00161，最大单步增量/初值约 0.00200），场探测面最大 |E|≈6.63×10⁷ V/m、最大 c|B|≈6.60×10⁷ V/m。`held-no-feedback-control` 下游最大轨道差约 2.54×10⁻¹⁵ µm；`held-no-initial-field` 的场能量从零增长至约 4.41×10⁻¹¹ J，因此仍可见零场启动效应。额外的强横向动量组（|u_x|、|u_y| 可达约 3）在启动前仍保持 x/y；无反作用的下游最大轨道差约 1.36×10⁻¹⁴ µm，场能量峰值/初值约 1.00161，最大单步增加/初值约 0.00161。步内启动的分段 Yee 电流测试中，冻结段 `J_x=J_y=0`，恢复动量后横向电流非零，离散连续性相对残差约 1.3×10⁻¹⁶。高 boost（γ≈929）的独立几何检查中，完整三维同步导致样例粒子的初始横向位移约 16.3 mm，横向保持则保留输入面 x/y；相同实验室 u_z、不同横向动量的冻结纵向速度一致，输入面时刻误差约 4.3×10⁻⁸ fs。以上只说明这个小型自由漂移数值回归没有检测到病态场能增长，**不能**推出真实高 boost、40 µm 束团的场源完整、无界面辐射或物理收敛；高 boost 的场能量本次未作生产尺度测试。仍须用实际 HDF5 扫描孔径、CPML、粒子权重、时间步长，并重点比较启动前后 E/B 与场能量的单步变化。冻结模式的能量账本标记 `closure_valid=false`，场能量比值只能作为诊断，不能当作能量守恒证明。

若现有运行同时设置 `mesh.ballistic_outside_transverse: true`，横向物理区外的冻结粒子在进入物理区前仍不向 Maxwell 网格沉积电流。初场对横向 CPML 粒子的可选一次性求解也不等于连续的框外场源；这部分误差仍需单独评估。
