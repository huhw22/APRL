# 场—粒子交换的时间、空间与轨迹配对诊断

本目录包含只读分析、独立诊断程序与默认关闭的完整程序观察钩子。
Cowan、Boris、电流沉积、CPML、供场的数值更新以及生产计算卡保持。
时间诊断结果在`results/particle_grid_exchange_diagnosis_20261031/r1/README.md`；
本轮空间／轨迹诊断与半隐式研究在
`results/particle_work_geometry_20261031/analysis-r1/README.md`。

`analyze_time_pairing.py`在计算节点汇总19099与19098已关闭的每rank严格账本，
直接量出P−M_mid=(P−M_old)+(M_old−M_mid)。P扣除原记录的边界、frozen、
carrier约束源功；old/mid使用同一份已沉积电流。结果是原轨迹的恒等分解，
不能当成修改时间取样后的因果预测。每rank完成标记、步号、有限值及
current_time_pairing恒等式均检查，原文件size/mtime保持。

`external_magnet_probe.cpp`只链接现有`libdirect_eb_core.a`，用同一boost、
原170周期波荡器与原dt，推进三个真实输入的单粒子，比较1/5/10/20个子步。
参照为静态lab磁场不做功；没有Maxwell更新、反馈、电流沉积或CPML。
它能排查抽样粒子的外场直接误差，不能冒充完整PIC或全束物理验收。
`run_external_probe.py`在计算节点编译独立可执行文件，核心库和头文件均冻结。

`present.py`仅用保存的小CSV/JSON作图和生成报告。主曲线保持既定峰值分母，
另有明确标注终态固定分母的比较；不同功项保留符号。两项相消变小不等于
逐项收敛。公式、代码时序、原数据和原型优先级在生成的报告中说明。

两个sbatch入口分别限10分钟/8CPU、5分钟/2CPU；都是独立诊断，
不启动全量PIC，不替换生产程序或默认入口。后续实际时间取样／轨迹规则
变更需要先确认，不能为了账本闭合自行缩放粒子能量。

## 默认关闭的空间与轨迹观察

`APRL_PARTICLE_WORK_GEOMETRY=/绝对目录`在所有rank一致设置后启用；
需要原能量账本、活动移动窗和自场反馈，无filter/checkpoint/retirement/
ballistic-outside/incident。目录须由运行者创建，既有输出拒绝覆盖。
每rank逐步CSV及complete标记；关闭时无逐步新增MPI交换，但有一次启动设置校验，
不声称严格零成本。原批处理Trace内有固定大小的诊断状态，不存全粒子轨迹历史。

启用时私有一层电场邻居延拓支持跨rank路径，缓冲不安装进物理网格。
独立二点Gauss求积匹配原电流沉积的原生空间权重，与原三线性点值取场分别量测。
CSV观察的是数值推送尝试，实际物理退出事后可能裁切；
`present_geometry.py`额外联用原v5移除计数筛选所有粒子释放后、首次移除前的完整段。
全段其它口径保留P−B和C−Mold，不能当作纯空间误差。

`geometry_checks.cpp`和`geometry_mpi_checks.cpp`分别验证3D原沉积伴随功与跨rank路径；
`prepare_geometry_cases.py`冻结短有磁／漂移的原程序、新程序off/on比较。
`complete_geometry_campaign.py`只读复验已关闭输出，并补启动负控和经典线性响应小原型。
`linear_response_probe.py`没有空间运动或Maxwell传播，不是全PIC半隐式实现。
研究方法和适配约束见[SEMI_IMPLICIT.md](SEMI_IMPLICIT.md)。
