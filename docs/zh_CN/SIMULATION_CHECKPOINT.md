# 模拟快照与续跑

模拟快照保存一个完整场步结束时的计算状态。后续运行可以从这个状态继续推进，
省去重新读取、变换整个初始束团以及再次求解初始自场的成本。快照采用同步写入：
所有 MPI rank 完成写入并发布文件后，原模拟继续运行。

快照时间是 **boost 计算系的共同时间切片**，不是整个束团处于相同 lab 时间的切片。
可以在任意完整 boost 场步保存：束头已经在下游、束尾仍在磁场中、部分粒子仍处于
held 状态，都不妨碍保存和精确续跑。无需等待整个束团离开波荡器，也不会因此自动
修改窗口、吸收层或物理停止位置。

## 保存指定步或指定时刻

在已有输入卡中增加顶层配置：

```yaml
units:
  length: m
  time: s

checkpoint:
  directory: /absolute/path/run-a/checkpoints
  steps: [1000, 2000]
  times: [1.0e-11]
  io_timeout_seconds: 3600
```

- `steps` 是从原模拟开始累计的完整场步数；`0` 表示初始化及初始诊断完成后的状态。
- `times` 是从原模拟开始累计的绝对 boost 时间，继承 `units.time`。目标落在两个场步
  之间时，在第一个达到或超过目标时间的完整步结束处保存；不插值，也不改变时间步。
- 步列表和时间列表分别排序、去重；同一步触发多个请求，只写一份快照。
- 可以只提供其中一个列表。省略整个 `checkpoint` 块，或两个列表均为空且不指定
  `restart_from`，不会创建快照文件。
- `io_timeout_seconds` 是一次快照读写的本地超时保护，默认 3600 秒，必须为正整数。
  它独立于原有 MPI 收尾保护，不修改正常场推进的超时策略。

文件名采用实际完成步，例如 `step-000000001000.h5`。日志中的
`[checkpoint] phase=write` 给出实际步数、实际 boost 时间、写入字节数、总等待时间、
有效 GiB/s 和文件路径。写入等待从推进耗时统计中单独扣除，不能把它计作每步推进成本。

未启用时，推进循环仅保留一个关闭的条件分支，不新增快照 I/O、MPI 集合通信或全场复制。
启用但尚未到达保存时刻时，只比较本地步号和时间。写入速度需以目标集群和文件系统上的
日志实测为准；文件系统负载会影响结果。

## 从快照精确续跑

复制原输入卡，为全部输出选择新的目录和文件路径，然后设置：

```yaml
output:
  overwrite: false
  manifest: /absolute/path/run-b/run-manifest.yaml

checkpoint:
  restart_from: /absolute/path/run-a/checkpoints/step-000000001000.h5
  directory: /absolute/path/run-b/checkpoints
  steps: [3000]
  allow_retirement_change: false
  reset_diagnostics: false
  io_timeout_seconds: 3600
```

`mesh.duration` 仍是从原模拟零时刻计算的**绝对 boost 时间上限**，不是续跑时长；它必须
大于快照保存时间。后续 `steps`、`times` 也保持原时间轴，恢复时跳过不晚于保存点的请求。
轨迹、探测器和能量账本的输出路径也必须改到新目录。续跑拒绝 `output.overwrite: true`。

默认要求物理模型及诊断定义与父输入卡一致。网格、boost、束流、磁场/源、边界、反馈、
初场配置和能量账本计算设置保持不变；可以调整输出路径、缓冲及压缩参数、运行报告设置、
`mesh.duration` 和停止条件。修改输入卡的写法不应被当作语义等价保证：兼容检查比较
规范化后的 YAML 定义，例如单位表达方式也应沿用原卡。

第一版要求相同 MPI rank 数和 z 分片，并检查本机算术表示兼容性，尤其是 `long double`
累积状态。快照不支持改 rank 数后的自动重分片。需要逐位比较时，应保持原工具链、计算
实现和 MPI/OpenMP 布局；换节点不等于改变 rank 数。

## 保存哪些状态

快照包含：

- 六个交错 E/B 场数组及 CPML 历史变量，保留原 leapfrog 时间交错；
- 全部本地粒子的顺序、ID、位置、动量、权重、held 状态等，以及 CPML 和退休电流载流子；
- 解析后的网格、分片、boost 原点、绝对时间、步数、粒子子步信息及累积边界计数；
- 保持原精度的能量账本累积和初始参考能量，以及轨迹和场探测面的下一采样时刻；
- ballistic 双平面验证尚未匹配的入口记录及已有累计统计；
- 原始输入卡、运行身份、父快照，以及诊断文件的历史前缀边界。

保存点位于完整场更新和本步诊断之后。下一步重新构造的临时电流/工作缓冲无需作为
历史状态保存。读取快照后不会重新做初始 Poisson 求解，也不会重放快照之前的粒子推进。

## 从同一状态比较未来的吸收方案

仅改变快照之后的粒子退休层，可以显式开启：

```yaml
checkpoint:
  restart_from: /absolute/path/run-a/checkpoints/step-000000001000.h5
  allow_retirement_change: true
  reset_diagnostics: false

particle_retirement:
  enabled: true
  entrance_z: 4.5
  length: 0.3
```

示例长度遵循本卡 `units.length`。该开关只允许修改 `particle_retirement`，不允许顺带
修改网格、磁场、粒子权重或积分方法。原有几何、磁场排斥和探测保护检查仍然执行。

新的退休分支必须从旧退休处理尚未发生的状态开始：不能已有退休入口计数或退休载流子；
若新方案启用退休，当前所有粒子还必须位于新入口上游。否则程序拒绝分支，要求选择更早
的快照。原因是后来修改吸收方法无法撤销已经写入场中的旧电流历史。

这一检查不要求所有粒子离开磁场，也不要求 held 全部解除。束头与束尾可以处于不同
阶段，只要新的退休处理尚未作用于任何粒子。自由传播长度、最终探测位置和何时停止仍由
输入卡的物理设计决定，快照功能不会自动延长它们。

## 新诊断段与近场参考记录

`reset_diagnostics: false` 恢复原采样相位，续跑输出从下一次原定采样开始，避免重复写入
快照时刻。能量账本继续沿用原初始能量和累计流出能量，不把重启点当成新的能量零点。

`reset_diagnostics: true` 允许改场探测面、粒子探测面和轨迹输出定义，并从当前时间创建
新的诊断段。当前在窗口内的场平面可从恢复时刻开始采样，轨迹也从该时刻建立新的节拍。
该选项不重置能量账本的物理累计，也不会生成快照之前的场或粒子穿越记录。

新增 ballistic reference 诊断尤其需要较早的快照：若有粒子已经穿过新参考入口，或此前
已有粒子丢失/退休，单靠当前状态不能补回真实入口事件，程序会拒绝这类缺历史的配置。
已有 reference 诊断采用精确恢复时，则保留其未匹配入口和累计验证状态。

ballistic reference 文件保存的是入口穿越记录；粒子穿越入口后，即使仍留在快照中，
续跑也不会再次生成那条入口记录。近场重构必须使用父段已有入口记录及续跑新增记录，
不能只读取新目录后把缺失的参考电荷解释为辐射信号。

## 历史前缀与后处理接续

写快照前，程序将轨迹、场/粒子探测器和能量账本的缓冲提交并刷新到文件，但不关闭这些
原输出。快照的 `metadata/diagnostic_prefixes` 为每个 MPI rank 保存一段 YAML，包含：

```yaml
step: 1000
time_box_s: 1.0e-11
parent_checkpoint: ""
outputs:
  - file: /absolute/path/run-a/detectors/field.h5
    group: /field_plane
    records: 101
    validation_records: 0
```

这些是各 rank 数据块中的元数据，不是把所有 YAML 混成一个文本；实际 HDF5 路径为
`/blocks/metadata/diagnostic_prefixes/data`，同组 `rank_counts` 给出各 rank 的字节数。
root 记录探测器和账本，每个 rank 记录自身轨迹。`records` 是
快照时已经提交的样本/行数，`validation_records` 是 reference 校验记录的独立前缀长度。

后处理接续规则：

1. 从续跑输出的 `restart_parent` 找到父快照。
2. 读取父快照给出的文件和前缀计数；只读取相应前缀。父运行后来可能继续写入，不能把
   父文件当前全部内容与新段拼接，否则会重复计入快照之后的区间。
3. 追加续跑段当前 `committed_samples`/`committed_records` 范围内的数据。
4. 如果父快照本身来自续跑，再递归接续其 `parent_checkpoint`。新诊断段只能与定义相容
   的历史段合并，不能把更改位置或节拍的文件当成同一个连续探测面。

续跑 HDF5 group 的 `restart_parent`、`restart_step`、`restart_time_box_s` 和
`output_scope` 明确标记这种分段关系。`complete=1` 表示本续跑段达到物理停止条件，
不表示该单个文件含有从原始 t=0 开始的全部历史。

仓库提供只读前缀清单工具 `postprocess/checkpoint_segments.py`。它只依赖 Python 的
`h5py` 与 `PyYAML`，无需 MPI 启动，也不依赖其他后处理程序：

```bash
python3 postprocess/checkpoint_segments.py \
  /absolute/path/run-b/checkpoints/step-000000003000.h5 \
  --verify-sources > checkpoint-prefixes.json
```

工具递归读取祖先快照，输出从最老一代到指定快照的 `generations` 与有序 `slices`。
每个 slice 明确给出 `generation`、`rank`、`file`、`group`、`records` 和
`validation_records`，其范围均为相应数据集的 `[0:count]`。`reset_diagnostics` 保留在
每代和每个 slice 中，工具不会把不同诊断定义的同名文件自动拼成一个探测面。

默认检查祖先环、重复文件/组、文件缺失、快照完成状态以及各 rank 的步数/时间一致性。
`--verify-sources` 进一步检查源 HDF5 的 `run_id`、续跑父身份、提交计数及数据集维度，
要求文件中已提交范围至少覆盖冻结前缀。检查只读取元数据、标量和维度，不加载、复制或
合并庞大的场数组。父模拟仍在写输出时，HDF5 可能暂时拒绝打开其文件；默认清单模式只
检查这些文件是否存在，详细源文件验证可在父输出正常关闭后执行。

输出只包含**截至指定快照**的历史，不自动搜索它之后产生的续跑文件。既有物理后处理器
仍不会自动拼接该 JSON；分析脚本可按 slice 清单分批读取，或显式构造合并数据。不要
删除/覆盖用于接续的原诊断文件和父快照。仅做动力学续跑时，大数组与累计状态已在快照
内部；历史诊断文件主要用于补全之前的输出。

对已有小型验收目录的只读自测可使用：

```bash
python3 postprocess/checkpoint_segments.py --self-test \
  /absolute/path/results/checkpoint_20261026/science-r1
```

自测读取步 7→17 的两代快照及跨入口/出口的 ballistic validation 文件，并用内存中的
异常元数据检查环、重复前缀、缺文件和超出已提交范围的拒绝行为；不会改写验收结果。

## 并行写入、完成标记与失败

快照使用并行 HDF5/MPI-I/O。每个 rank 直接写入自己拥有的数据区，不将全场或全粒子
集中到 rank 0。大数组无压缩写入，使用有界传输批次；小型标量和状态元数据另存。

写入期间文件名带 `.partial` 后缀。所有 rank 成功关闭文件后才发布最终 `.h5` 文件；
恢复只接受通过格式、布局及完成状态检查的最终文件。已有最终文件或 `.partial` 文件
都不会被覆盖，即使普通输出允许覆盖也是如此。失败的 `.partial` 留作诊断，不能当成
可恢复快照使用。

I/O 期间沿用本地、无需 MPI 配合的 watchdog 思路，超时或任一 rank 报错将使本次作业
失败，防止某个节点退出后其他节点永久停在通信阶段。普通的 MPI 收尾保护仍然保留。
集群吞吐和数值等价性验收结果应以具体测试目录的记录为准；本说明不替代实测。
