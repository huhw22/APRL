# 安装、卸载与依赖维护

本项目提供四个顶层维护脚本。它们只处理依赖、构建目录和已安装副本，不提供删除源码树的操作。

## 支持的编译环境

核心程序的最低接口要求为：

| 组件 | 要求 | 说明 |
|---|---|---|
| 操作系统 | Linux/Unix 风格环境，Bash | 自动包安装器目前面向 Linux |
| CMake | 3.16 或更新 | 安装和依赖探针使用相同下限 |
| C++ | 完整 C++11 编译器 | 必须与 MPI wrapper 所用编译器 ABI 一致 |
| MPI | 提供 C++ 编译 wrapper 和运行器 | 不限定 Open MPI/MPICH；编译和运行必须使用同一家族 |
| HDF5 C | 1.10 或更新 | 生产 MPI 输入要求 parallel HDF5 |
| yaml-cpp | 0.6 或更新 | 主程序和多个独立工具需要 |
| FFTW3 | double precision；完整安装需要 `fftw3_threads` | 仅 `full` 后处理配置需要 |
| pkg-config | 可解析 FFTW3 | 仅 `full` 配置需要 |

本机端到端安装测试使用 GCC 13.3、Open MPI 4.1.6、parallel HDF5 1.10.10、yaml-cpp 0.8.0 和 FFTW 3.3.10。这是已验证组合，不表示程序只能使用这些精确版本。

原生 Elegant SDDS5 转换器是可选组件。它要求已经编译好的官方 SDDS 源码树，并同时提供 `SDDS1`、`rpnlib`、`mdbmth`、`mdblib`、zlib、liblzma 和 GSL。脚本不会用一个不完整的自写 SDDS 解析器代替这些库，也不会在后台下载 SDDS 源码。

### 最容易出现的交叉依赖错误

MPI C++ wrapper、MPI launcher 和 parallel HDF5 必须来自同一套 MPI/编译器环境。例如，用 Open MPI 的 `mpic++` 链接由 MPICH 编译的 HDF5，即使 CMake 能找到文件，也不是可用环境。集群上还应避免同时启用系统库、Conda MPI 和管理员 module 中的 HDF5。

建议从干净的登录环境开始，只加载一套编译器/MPI/HDF5 module，再运行依赖体检。`CMAKE_PREFIX_PATH`、`HDF5_ROOT`、`PKG_CONFIG_PATH` 和 `SDDS_ROOT` 可以用于明确选择管理员提供的非系统安装。若需要额外 CMake 参数，安装脚本支持重复的 `--cmake-arg`。

## 依赖体检

完整检查：

```bash
./check_dependencies.sh --profile full --report dependency-report.log
```

若已有官方 SDDS 构建：

```bash
./check_dependencies.sh --profile full \
  --sdds-root /path/to/built/SDDS \
  --require-sdds \
  --report dependency-report.log
```

检查器执行的不是简单的命令名查找，而是：

1. 记录操作系统、工具版本、MPI wrapper 和 HDF5 配置信息；
2. 让 CMake 重新解析 MPI、HDF5、yaml-cpp、FFTW 和可选 SDDS；
3. 编译、链接一个混合依赖探针；
4. 对比 HDF5 编译头文件与运行时库版本；
5. 执行单 rank 测试；
6. 执行两 rank collective parallel-HDF5 写入测试；
7. 记录最终动态库解析路径，便于发现系统库和 module/Conda 库混用。

报告开头始终包含管理员可读的依赖清单。即使 CMake、网络或某个库不可用，已有日志仍可直接交给集群管理员。

某些集群禁止在登录节点直接启动 MPI。此时可以先执行：

```bash
./check_dependencies.sh --profile full --skip-mpi-run \
  --report dependency-report-login-node.log
```

然后在 `salloc`/`sbatch` 分配到的节点内重新运行不带 `--skip-mpi-run` 的命令。跳过测试只表示 MPI 运行阶段未验证，不应当被解释成通过。
若站点要求 `srun`，可使用 `--mpi-launcher srun`；指定的运行器必须支持通用的 `-n 2` rank 参数。

`--allow-serial-hdf5` 仅用于单 rank 本地调试。默认配置会拒绝 serial HDF5，防止程序成功编译后才在超算多 rank 读粒子时失败。

## 自动安装系统依赖

在可联网、允许包管理器访问的软件机上，以管理员身份执行：

```bash
sudo ./install_dependencies.sh --profile full \
  --report dependency-install.log
```

工具支持 `apt-get`、`dnf`、`yum`、`zypper` 和 `pacman`，安装编译器、CMake、MPI、parallel HDF5、yaml-cpp、FFTW 及 SDDS 所需的系统支持库。包名随发行版和软件源可能变化；若安装失败，完整的包管理器输出仍保留在日志中。

这个脚本有意要求 `root`，不会在内部自行调用 `sudo`，也不会修改源码。它不负责安装官方 SDDS。安装结束后会以 root-safe 模式重新编译和链接依赖探针，但跳过两 rank 启动；普通用户仍应在实际作业环境中完成 MPI 运行检查。

离线超算通常不应运行自动包安装器。直接运行 `check_dependencies.sh` 生成需求日志，并让管理员提供一致的 module。若站点有离线软件源，可以由管理员使用 `--skip-refresh` 调用安装器。

## 安装程序

默认安装完整配置，包括主程序、通用小工具和五个后处理器：

```bash
./install.sh --prefix /path/you/can/write
```

若需要系统安装：

```bash
sudo ./install.sh --prefix /opt/aprl
```

如果需要原生 SDDS 转换器：

```bash
./install.sh --prefix /path/you/can/write \
  --sdds-root /path/to/built/SDDS \
  --require-sdds
```

安装顺序固定为：依赖体检、全新构建目录、统一 Release 编译、必需回归测试、安装、写入回执、删除编译结果。`full` 配置把五个后处理器纳入同一个构建和端到端测试图；`core` 配置显式关闭它们。任何依赖、编译或测试失败都会在安装前终止，并保留失败构建供检查；只有成功安装后才自动清除构建目录。

常用选项：

- `--profile core`：只构建主程序和核心小工具，不要求 FFTW；
- `--keep-build`：成功后保留编译目录；
- `--skip-tests`：跳过回归，仅适合已经在完全相同环境验证过的重复打包；
- `--skip-mpi-run`：登录节点不允许 MPI 启动时使用；
- `--mpi-cxx /path/to/mpic++`：明确 MPI wrapper；
- `--mpi-launcher srun`：明确两 rank 依赖检查使用的运行器；
- `--generator NAME`：明确 CMake generator；
- `--cmake-arg ARG`：传递站点 toolchain 或依赖路径，可重复；
- `--allow-serial-hdf5`：只构建单 rank 本地版本，并跳过 parallel-HDF5 专项测试。

默认前缀是 `/usr/local`，默认一次性构建目录是源码下的 `build-install`。安装完成后，二进制位于 `<prefix>/bin`，文档、示例和绝对路径安装回执位于 `<prefix>/share/aprl`。其中主程序为 `<prefix>/bin/aprl`。

重复安装到同一前缀时，脚本会先保存并验证旧回执；新安装及新回执成功后，才删除旧版本中已经不存在的文件。这样更新不会因某个已移除的工具留下不可追踪的旧二进制。

## 卸载

先检查目标：

```bash
./uninstall.sh --prefix /path/you/can/write --dry-run
```

确认后卸载：

```bash
./uninstall.sh --prefix /path/you/can/write
```

卸载器只读取 `install.sh` 创建的绝对路径回执，并拒绝任何不在指定前缀内的条目。它删除安装副本和默认的 `build-install` 编译目录，但不删除输出数据、依赖库、系统包或源码。没有合法回执时不会根据文件名猜测并删除内容。

系统前缀下卸载需要相同权限：

```bash
sudo ./uninstall.sh --prefix /opt/aprl
```

自动依赖安装器安装的是系统共享软件包，因此普通程序卸载不会反向删除它们；这些库可能仍被其他程序使用。
