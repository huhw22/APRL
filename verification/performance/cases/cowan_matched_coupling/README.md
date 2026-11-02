# Cowan 配套自场取样的独立候选

目的：在保留 Cowan 场推进、原粒子动量推进、守恒电流沉积、参考供场、释放事件和 CPML 的前提下，使粒子电功与 Cowan 的配套能量、纵向磁冲量与配套场动量一致。实现是实验覆盖层，不替换生产程序和默认配置。

`metric_gather.h` 给出局部取场算子；`native_metric_halo.h` 提供其必需的纵向邻层；`prepare.py` 从冻结的 `coupled_sbp_r2` 源码/库建立独立程序，仅重编 `source_coupling.cpp.o`。整个 Maxwell 网格不作新增迭代，也不预测未来粒子状态。MPI 消息数量、标签和迭代集合调用保持，邻层载荷增加；不能据此声称零成本。

数学关系及限制见 [METHOD.md](METHOD.md)。程序验证和物理对照分开：

* `checks.cpp`：真实沉积器、独立嵌套九点权重的功和磁冲量恒等式，以及省略权重的负控。1/2/4 rank 覆盖双向跨分片轨道。这是实现验证。
* `periodic.cpp`：原 Cowan/Ampere 核与真实带电粒子、原沉积器构成闭合周期盒。正负粒子初始重合，总电荷为零；有非零速度、电流和场粒子能量交换。保留原生 E/B 能量，并独立测量配套能量、场动量恒等式、电/磁粒子冲量余差。轴向波和斜向波分开；同一网格改变 quiet-start 的粒子采样数，**不是网格收敛测试**。
* `prepare_physics.py`：逐项复用 `source_coupling_exit_20261031/physics-r1` 的冷均匀球、会聚均匀球、冷高斯输入和理论到达事件；仅重定向输出并选择新程序。原严格账本、粒子轨迹及三个理论诊断面保持。每例 2700 秒硬限，单节点 16 MPI×8 线程，不指定节点。旧理论、释放相空间和电荷不作重新拟合。

当前证据入口为工作区 `results/cowan_matched_coupling_20261102/README.md`。粗加载的周期控制已暴露明显电力冲量余差，不能隐藏该失败例，也不能把更密加载下的小余差或数学能量闭合直接称为普遍物理验收。全程序移动边界、CPML 和实际束团的结论须由对应物理对照给出；原生场输出和真实辐射精度仍需单独核对。

构建和短验证：

```sh
python prepare.py --workspace /absolute/workspace --output /absolute/fresh-result --build /absolute/fresh-build
sbatch /absolute/fresh-result/run.sbatch
```

空间电荷重放：

```sh
python prepare_physics.py --workspace /absolute/workspace --build-record /absolute/closed-build-record --output /absolute/fresh-physics-result
sbatch /absolute/fresh-physics-result/scripts/physics.sbatch /absolute/fresh-physics-result uniform
```

同理分别提交 `converging` 和 `gaussian`。确认空闲资源后提交；不得改运行中冻结卡，不读未闭合 H5。结束时读取各组 `aprl.log`、组级 `*_COMPLETE.json`、`analysis/RESULT.json` 和原始账本；执行完整性、动力学误差、数值收支误差及辐射适用性分别报告。

`prepare_mpi.py` 单独冻结双节点四rank的真实邻层检查。`prepare_flow.py` / `run_flow.py` 重放原257步有磁/无磁小卡，总16 MPI×8线程分布到两节点；其旧参考场账本明确标记未采用新配套度量。`review.py` 核对冻结文件、原输入与理论身份、其它内核源码及原生产程序身份，汇总已闭合JSON和逐步CSV并生成相对图；不读取运行中H5。

`prepare_post.py` 在独立输出树对已闭合双节点证据重新汇总，不运行 PIC。首次19245因旧读取器省略新功差而FAILED；原记录和门限保留，19246仅修读取器加入独立功差项。原生物理收支和差额不归零。修复原因与公式见 METHOD.md。
