# 修正后全量漂移的只读统计

面向19098：完整41732775个输入、两面原始粒子、两面标量四动量输运、
boost严格账本，以及0.5 m原生lab E/B。结果在工作区
`results/corrected_full_comparison_20261030/free_drift_review_19098_r1/README.md`。

原模拟步正常结束，旧收尾检查器误把继承自SDDS的稀疏source_id当成1..N。
本工具只接受明确这一项报错且场已完整关闭的情形；按particle_id索引原输入行，
检查每面唯一完整覆盖，并逐行核对真正的source_id。原FAILED状态和旧完成记录保留。
通用check_completion.py已修正源模板，运行中的冻结脚本不改写。

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
