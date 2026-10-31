# 修正后全量漂移和波荡器的只读统计

面向19098：完整41732775个输入、两面原始粒子、两面标量四动量输运、
boost严格账本，以及0.5 m原生lab E/B。结果在工作区
`results/corrected_full_comparison_20261030/free_drift_review_19098_r1/README.md`。

原模拟步正常结束，旧收尾检查器误把继承自SDDS的稀疏source_id当成1..N。
本工具只接受明确这一项报错且场已完整关闭的情形；按particle_id索引原输入行，
检查每面唯一完整覆盖，并逐行核对真正的source_id。原FAILED状态和旧完成记录保留。
通用check_completion.py已修正源模板。19099实际收尾检查器的已授权替换另有
postcheck-amendments/source-id-rule.json记录；原冻结清单和原脚本仍保留，
本分析使用能识别这一份修正的validate_frozen，不扩大物理依赖的白名单。

`analyze.py`仅在计算节点运行；冻结同目录report_discrete.py/review_discrete.py、
relative_ledger.py、原FFT核analyze_spectra.py和独立小谱读取器independent_review.py。
逐粒子gamma差用差平方配对的稳定形式，电荷作为权重，std为总体标准差。
两面的标量端点能量用高低位和70位十进制数相减，保留两种诊断口径的差。
严格账本维持原固定D定义。FFT直接复用已验核，不扣场、不改时间层或Hann定义。

`run.sbatch`限制25分钟，一可用独占节点、16 CPU任务、128 GiB；不是新PIC。
`render.py`仅使用小CSV/JSON/NPZ画科研静态图与Markdown；两份旧漂移谱复用缓存。
`--reuse-figures`只补充文字和指标，保持已核验图像不改写。

局部谱底用90.5–91.3与91.7–92.5 eV的原生频格做线性拟合，
对91.3–91.7 eV原生频格求有符号超额，不截负、不改变数组；
这不是近场分离或纯辐射拟合。boost离散闭合也不构成完整lab世界管或动力学精度验收。

## 19099的部分出口覆盖

`analyze_undulator.py`针对已关闭的19099：完整原生场和32-rank严格账本，
入口41732775条、出口41732699条。它只接受这次明确的出口覆盖报错，
不把任意FAILED任务当成功；原检查器、FAILED和覆盖失败标记均保留。
共享scan_particles默认仍要求全覆盖；仅显式allow_exit_subset分支对真实共同ID
计算两面统计，入口同步排除缺失出口ID，另存全部入口统计及76个缺失ID清单。
标量端点也减去这76个ID的入口四动量，不给未知出口补零或外推。

`check_subset.py`用重排ID、任意source_id及已知能损验证共同cohort，
同时检查严格模式拒绝缺失、重复和source错配。`run_undulator.sbatch`仍是
25分钟只读后处理；大H5只在计算分配内扫描，不运行新PIC。

`render_undulator.py`复用旧有磁18932/18992和当前漂移19098的小谱缓存，
保持各自完整时间记录与原生频格。既画既定峰值分母，也把同一交换差画成
相对终态净boost粒子变化，防止外磁场的大幅中途能量摆动掩盖小净交换误差。
结果在`results/corrected_full_comparison_20261030/undulator_review_19099_r1/README.md`；
完整分析不等于物理能量交换或纯辐射分离通过。
