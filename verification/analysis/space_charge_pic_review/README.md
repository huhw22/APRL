# 球形空间电荷 PIC 结果审阅

`review.py` 只读闭合的三组实际粒子面，独立使用 long-double 线性洛伦兹
变换复核原 light-cone 分析；区分逐粒子自身入口事件、实际 PIC 事件时刻
与参考轨迹自身过同一 lab 面的时刻。计算实验室能量偏差的静止系动能／
纵向动量分解，保留原始数值；动量重置列仅为诊断，不作为纠正后的结果。
读取器检查18面完整性、ID、电荷权重、有限值，并核验数据前后哈希。

```bash
python3 verification/analysis/space_charge_pic_review/review.py \
  --workspace /public/home/wfollowm/APRL_Workspace --output /absolute/new/review
python3 verification/analysis/space_charge_pic_review/render.py /absolute/new/review
```

输出逐面 CSV、汇总 JSON、输入 SHA256、验证记录、两张 PNG/PDF。
`SUMMARY.json` 中的 `ledger` 保持 boost 系口径，`*_rest_*` 为束团初始平均
静止惯性系；两者不同。相对账本分母沿用每组固定 max|delta K_box|。
三组相同电荷/初始rms但输入并非同一个同时刻释放：过面世界线已包含
进入 PIC 之前的参考演化，不把本轮当成零内场启动或完整冷起点验收。

基准参考为球对称交叉前的相对论壳层 ODE。初始球均匀或高斯，
后者并不被约束始终保持高斯。没有按结果调节理论幅度、拟合空间电荷
强度或修改 PIC 推进；运行完整、账本闭合与物理符合分别报告。
