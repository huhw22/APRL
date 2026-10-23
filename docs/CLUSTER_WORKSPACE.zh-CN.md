# APRL 集群统一工作区

## 打开哪个目录

在远程 Codex 项目中选择整个 `/public/home/wfollowm/APRL_Workspace`。源码、实验脚本、输入链接、编译产物和结果都从这个根目录管理；不要只打开 `source/aprl`。

```text
APRL_Workspace/
├── AGENTS.md                     # 新会话必须遵守的项目约束
├── HANDOFF.md                    # 本次交付版本、状态和下一步
├── WORKSPACE_LAYOUT.md           # 本说明的工作区副本
├── cluster.env                  # 服务器路径和运行环境；source 后使用
├── source/
│   ├── aprl/                     # 候选源码，独立 Git checkout
│   └── aprl-baseline/            # 7a8ff35 基线，独立 Git checkout
├── build/
│   ├── candidate/               # 候选构建产物
│   ├── baseline/                # 基线构建产物
│   └── verification/            # 如需完整 CTest，另建此目录
├── tests/
│   ├── performance/             # 小型对照提交入口及 Slurm 日志
│   ├── elegant_4m/              # ELEGANT 实验脚本及独立 runs/
│   └── sase_4m/                 # 原 4 m SASE 实验脚本
├── inputs/
│   └── elegant-particles.h5     # 指向服务器已有输入的只读用途链接
├── results/                     # 本工作区新产生的结果
└── logs/                        # 构建和其他检查记录
```

`source/aprl/verification` 保留在源码仓库中，包含版本化回归测试和实验模板。实际提交时使用 `tests/` 的工作副本，编译使用 `build/`，输出使用 `results/` 或对应工况自己的 `runs/`。这样实验文件不会混入源码修改。

## 旧目录如何处理

旧源码 `/public/home/wfollowm/aprl_dev`、旧实验 `/public/home/wfollowm/APRL_Project`、旧安装 `/public/home/wfollowm/aprl_dev/aprl-runtime-lab` 保留。本次交付不对这些目录执行移动、清理或覆盖。

其中旧源码现场 HEAD 为 `b37603f5e56215d6227882f2b8b5aa2b44d48de1`，与本次性能对照基线 `7a8ff35f5784df8868f24efe56f24132e4c5e678` 不同。不能把旧安装二进制直接标成 `7a8ff35`；请从 `source/aprl-baseline` 用同一服务器工具链重新构建。

输入链接的目标是 `/public/home/wfollowm/APRL_Project/APRL_Elegant_4m_2nodes/particles.h5`。交付时该文件约 2.67 GB，未复制到 Git，也未重复占用一份数据空间。链接不代表数据已迁移；不得删除旧目标或通过链接修改输入。每次正式运行记录该输入的路径、大小、修改时间和必要的哈希。若以后需要完全独立的数据区，应先核对可用空间，再复制并校验，最后切换链接；不要在运行中更换输入。

Git 远程只保护已经提交的源码和脚本。未提交修改、粒子数据、历史结果、Python 环境、编译产物不因存在 Git 远程而自动有备份。整理历史实验时先核对正在运行的作业及其绝对路径。

## 如何维护脚本

- 主源码及模板的 Git 远程：`git@git.tsinghua.edu.cn:huhw22/aprl_dev.git`。
- 本次候选分支：`wfollowm_dev/compute-performance-20261023`。
- ELEGANT 模板：`source/aprl/verification/performance/cases/elegant_4m`。
- SASE 模板：`source/aprl/verification/performance/cases/sase_4m`。
- 在 `tests/` 试运行时只写该工况和结果目录；验证后的脚本改进应回写相应版本化模板并提交 Git。
- 更新模板到 `tests/` 时只更新脚本/配置，先对比本地变更；不整目录覆盖 `particles.h5`、`runs/`、`work/`、`.venv/` 或日志。
- 候选程序、基线程序和输入生成器一律使用明确的绝对路径，不依赖旧安装在 `PATH` 中的优先级。

工作区 `cluster.env` 覆盖原实验脚本的旧安装默认值。ELEGANT 仍把运行目录建在 `tests/elegant_4m/runs/`；SASE 的 `CAMPAIGN_WORK_DIR` 指向 `results/sase_4m`，其日志和 `.venv` 可以留在测试区。这些是路径安排，不改变物理输入。

## 远程 Codex 连接说明

2026-10-23 已修复桌面连接的 `socket hang up`：服务器启用了 `KillUserProcesses`，普通 SSH 会话退出会清理该会话启动的 Codex。现在使用账号的 `codex-ssh-app-server.service` 用户服务，服务文件在 `~/.config/systemd/user/`，已设置登录时启动；没有更改系统级清理策略。

可只读检查：

```bash
systemctl --user status codex-ssh-app-server.service --no-pager
codex app-server daemon version
```

新会话运行在哪台机器与聊天是否相同是两个概念。先核对 `pwd`、主机和工作区内 Git HEAD，再运行编译或提交。不要把本机 `/data/L_FEL_program` 的绝对路径用于集群。
