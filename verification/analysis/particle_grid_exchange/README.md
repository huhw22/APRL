# 场—粒子交换的时间配对诊断

本目录只新增只读分析和链接已有核心库的独立诊断程序，不修改Cowan、Boris、
电流沉积、CPML、供场或生产计算卡。实际结果位于工作区
`results/particle_grid_exchange_diagnosis_20261031/r1/README.md`。

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
