离散能量账本（2026-10-30 验收收尾）

用途
  检查实际三维 Cowan/粒子/移窗步骤是否符合其离散能量收支。
  新账本与原 v5 HDF5 连续公式监测并存。推进、Boris、沉积、CPML、
  参考供场和原 MPI 更新顺序保持；原 H5 的字段和 closure_valid 不改。
  不是把新残差当成物理精度证明，也不把数值项当辐射/真实外功。

程序
  /public/home/wfollowm/APRL_Workspace/build/discrete-energy-ledger-20261029/aprl
  SHA256 c48c376929fbdabde4f9cbf08145aa11ea83f86e11839aec24e5f26f3601977b
  这是独立已验程序，默认 cluster.env/旧生产程序未替换。

启用
  在现有计算卡 energy_ledger 下加入：

    energy_ledger:
      enabled: true
      discrete_operator_audit: true
      directory: /你的结果目录/ledger

  其余原配置保持。默认 discrete_operator_audit=false；不开启时没有
  新场数组、场遍历、CSV或额外 MPI 面交换。关闭的严格耗时差没有单独
  做新旧二进制配对，不能声称数学意义上的零成本。

当前支持范围
  moving window + 原账本启用 + 未滤波 current + 原横向CPML物理孔径。
  checkpoint/restart、incident waves、ballistic_outside_transverse 组合拒绝；
  不改变这些功能的默认旧行为。物理粒子、held/快速参考等推进原样。
  独立缓冲仅覆盖活动窗 E/B，不恢复原来的长窗口。

输出
  ledger/discrete-field-rank-00000.csv 等：每rank、每完整场步的独立数值项。
  ledger/discrete-field-rank-00000.complete 等：完整标记和场步数。
  原 ledger/energy-ledger.h5 仍为 v5，原探测面/轨迹输出保持。
  本轮默认不覆盖已经存在的同名新CSV；复算应使用新结果目录。

完成后汇总（普通生产目录也能用，无需off控制工况）
  source /public/home/wfollowm/APRL_Workspace/cluster.env
  "$PYTHON_BIN" \
    "$APRL_SOURCE_DIR/verification/performance/cases/discrete_energy_ledger/report.py" \
    --ledger-directory /你的结果目录/ledger \
    --output /一个单独的空分析目录

  若旧H5不是energy-ledger.h5，加 --physical-ledger /路径/自定义文件.h5。
  汇总核对H5的mpi_size、完整性、全rank/全步覆盖、算子逐点预测误差、
  数值恒等式；支持旧H5稀疏采样步长。输出summary.json、discrete_energy.csv、
  discrete_energy.png/pdf。大文件后处理仍应放计算节点；本轮只读小工况。

主要字段的含义
  native_*                 原生Yee样本平方的场能，控制体积剔除横向CPML。
  discrete_curl_exchange    配对curl得到的边界/接口交换，不是远场辐射能流。
  applied_current_mid      只计电场更新实际采用的J，与前后E平均配对。
  unapplied_current_mid    最外层保持的横向E面上未进入Ampere更新的J功，单列。
  time_electric/magnetic   实际场增量的平方项；不是粒子或外场真实做功。
  current_time_pairing     applied_current_mid - applied_current_old。
  cowan_operator           Cowan与普通Yee Faraday算子功的差，显式保留。
  magnetic/electric_cpml   对所统计物理孔径的直接变化，非整个CPML区域耗能。
  inlet_overwrite/remap    原生能量口径下直接实测的覆盖、移窗变化。
  field_operator_residual  以上项目的独立核对差，需与预测误差一起检查。
  max_*prediction_error    按输入字段独立重算离散增量，与实际E/B增量比较。

新账本对上的含义
  表明当前实现与其离散算子/边界操作相符，并量化粒子-网格交换差。
  Cowan、时间层、空间定义和通量定义贡献没有被抹去；有些是表示方式，
  有些反映数值算子近似。物理精度仍需独立收敛/物理对照，不能据此认定
  所有近场/辐射分离、实验能损或lab系世界管守恒已验收。

本轮验证
  19027/node3：77 required；1/4 MPI三维制造场、有J、已知平移脉冲、
  故意改E的负控；原无磁257步同16MPI×8的off/on完整对照。
  19028/node3,5：4MPI跨节点同模块/负控；16MPI×8的off/on完整对照。
  所有原场/粒子面/轨迹逐位，v5科学账本字段逐位（计时字段分开）。
  新增B面交换只在启用时发生；仍MPI_FUNNELED，无新增停止barrier；
  原MPI-free收尾watchdog和600/60秒设置保持。
  reader另验7步采样与逐条数据一致，拒绝缺rank/未完成/损坏数值项/错误系。

结果（boost系）
  原残差                 -3.02832787249e-17 J
  新总核对差             -4.21410767111e-20 J（abs净粒子变化的0.0445737%）
  最大累计场算子残差     约1.13e-31 J
  零粒子/载流子/J尾段    旧-2.150087e-18 J → 新约-7e-33 J
  新总差主要是独立粒子-网格交换差；不是设定的通过百分比门槛。

原残差与新总差之间的独立数值/定义分解
  空间定义     -1.29047548308e-18 J
  通量定义     -2.94267197596e-17 J
  时间/电流配对+1.97982174160e-17 J
  Cowan算子    -1.93221598215e-17 J
  物理孔径CPML约舍入级。所有项独立计算，不从残差反推。

成本边界
  单节点40–240步：0.13650→0.25695 s/步，+88.24%；
  双节点同布局：0.13387→0.24689 s/步，+84.43%。
  每种分配各一对、无profile、完整小卡输出，不外推生产/严格关闭成本。
  最大观测总RSS分别1882.50→2245.71 MiB、1875.33→2252.36 MiB；
  这是日志采样到的并发RSS，不把各rank独立峰值相加当瞬时峰值。
  严格审计默认关闭，适合正确性验收和排错；未自动切换生产配置。

完整结果
  results/discrete_energy_ledger_20261029/
  campaign-r2/、mpi-r1/保存闭合原始证据；report-single-final/、
  report-mpi-final/保存通用汇总。campaign-r1仅编译告警清理前的未运行准备。
  原始运行后仅把review的“same-node”文字更正为“same-allocation”，数值
  逐字段未变，原报告保留为DISCRETE_REVIEW.original-run.json。
