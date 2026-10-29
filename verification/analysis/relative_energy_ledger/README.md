# 相对能量账本图与数据

主图与主 CSV 显示 `100 * signed_energy / D`，其中含粒子工况的
`D=max(abs(delta_particle_kinetic_box_J))`，取整个已完成记录的固定值。
ΔK是同一粒子集合含累计移除事件的稳定动能差；各条曲线共用分母。
100%对应本次动能变化的峰值尺度。该归一化不改变物理判据。

粒子动能变化严格为零时，依次使用初始、峰值原生场能；全零工况明确拒绝
定义相对量。图保留符号；不以每条曲线自身幅度缩放，也不逐时刻除以可能
过零的ΔK。`normalization.json`保存所用分母、时段和定义。

`energy_budget.csv`是百分数表；`energy_budget_absolute.csv`保留对应J值。
图的三个面板分别展示：粒子动能变化与独立场／边界／源项合计；总核对差
与粒子—网格交换差；场算子核对差。数值通道继续单列，不解释成辐射。

独立通道合计直接从场、边界、源功和数值项求出，不从残差反推。
当前期望的核对式为 `ΔK - channel_sum = total_residual`。

```bash
python3 verification/analysis/relative_energy_ledger/render.py --plan /absolute/plan.json
python3 -m unittest discover -s verification/analysis/relative_energy_ledger -p 'test_relative.py'
```

计划记录每组的关闭标记、标记中的成功字段、ledger路径和显示名称。
只有存在成功关闭记录的工况才会打开 H5；随后再次检查 H5 和逐rank CSV
完整性。没有关闭记录的运行中工况只列等待，不打开其H5。
复制到独立结果目录时，将既有`discrete_energy_ledger/report.py`另存为
`report_discrete.py`，并复制其`review_discrete.py`依赖。

当前19031使用的冻结脚本不改动。独立的依赖作业在19031结束后补绘其余
已完成工况；若原物理作业失败，只处理已关闭且完整的工况。
未来小卡模板已接入本模块，图和主数据默认使用相对量。
