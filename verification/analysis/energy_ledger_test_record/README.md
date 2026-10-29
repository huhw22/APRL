# 能量账本成功测试整理

用户交付为工作区 `results/energy_ledger_test_record_20261030/README.md`。
展示正文介绍已知真空脉冲的解析对照和含粒子自由漂移的能量收支。
自由漂移结果的结论限于离散能量收支一致性，不能扩展成束团动力学精度证明。
三维制造场、扰动识别、开关逐位比较、软件回归、读取器和成本检查
已移入独立的 `PROGRAM_VALIDATION.md`，不作为物理正确性的依据。

`assemble.py` 是一次性整理脚本：同文件系统移动已完成的
`results/discrete_energy_ledger_20261029/`，保留原路径的相对符号链接，
复制小脚本快照，并从已完成的 CSV/JSON 生成图片及 Markdown。
它不启动计算、不读取大型场数据重新分析、不改写原始证据；
目标目录已存在或原路径已经是链接时会拒绝重复执行。
其它诊断阶段的目录不在本次迁移范围内。

整理前应确认账号没有运行中的相关作业；`--confirmed-idle` 记录此检查，
不替代实际查询。首次执行示例（已完成的整理不要再次执行）：

```bash
python3 source/aprl/verification/analysis/energy_ledger_test_record/assemble.py \
  --workspace /public/home/wfollowm/APRL_Workspace \
  --target /public/home/wfollowm/APRL_Workspace/results/energy_ledger_test_record_20261030 \
  --confirmed-idle
```

只验证现有交付，不移动文件、不运行仿真：

```bash
python3 source/aprl/verification/analysis/energy_ledger_test_record/validate.py \
  /public/home/wfollowm/APRL_Workspace/results/energy_ledger_test_record_20261030
```

校验包括文件设备号／inode、非目录大小和 mtime、关键小文件 SHA256、
链接、脚本快照及程序哈希。大型 H5 使用同文件系统原位重命名保留身份，
没有声称重新扫描全部 H5 内容。完整原始验收记录保持原样。

现有目录只刷新两份文档及交付哈希时，可导入 `assemble.py` 的
`render_documents(evidence_root, bundle_root)`；不要重新执行搬迁入口。

依赖：Python、NumPy、Matplotlib；工作区 `physics` 环境已提供。
图的数值来自原始模块 summary 和逐步 CSV，绘图选用当前验收字段。
程序、测试卡和仿真核心在本轮没有修改。
