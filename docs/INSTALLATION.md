# 安装、卸载与依赖维护

本项目提供依赖引导、体检、系统包安装、程序安装和卸载脚本。它们只处理依赖、构建目录和已安装副本，不提供删除源码树的操作。

## 支持的编译环境

核心程序的最低接口要求为：

| 组件 | 要求 | 说明 |
|---|---|---|
| 操作系统 | Linux，Bash | 自动包安装器目前面向 Linux |
| CMake | APRL 本体 3.16；依赖引导 3.18 | HDF5 1.14.6 的源码构建需要 3.18；联网脚本会在版本不足时准备私有 CMake |
| C/C++ | C 编译器和完整 C++11 编译器 | C 用于 HDF5 配置探针；C++ 必须与 MPI wrapper 所用编译器 ABI 一致 |
| MPI | 提供 C++ 编译 wrapper 和运行器 | 不限定 Open MPI/MPICH；编译和运行必须使用同一家族 |
| HDF5 C | 1.10 或更新 | 生产 MPI 输入要求 parallel HDF5 |
| yaml-cpp | 0.6 或更新 | 主程序和多个独立工具需要 |
| FFTW3 | double precision；完整安装需要 `fftw3_threads` | 仅 `full` 后处理配置需要 |
| pkg-config | 可解析 FFTW3 | 仅 `full` 配置需要 |

本机端到端安装测试使用 GCC 13.3、Open MPI 4.1.6、parallel HDF5 1.10.10、yaml-cpp 0.8.0 和 FFTW 3.3.10。这是已验证组合，不表示程序只能使用这些精确版本。

原生 Elegant SDDS5 转换器是可选组件。它要求已经编译好的官方 SDDS 源码树，并同时提供 `SDDS1`、`rpnlib`、`mdbmth`、`mdblib`、zlib、liblzma 和 GSL。脚本不会用一个不完整的自写 SDDS 解析器代替这些库，也不会在后台下载 SDDS 源码。

### 最容易出现的交叉依赖错误

MPI C++ wrapper、MPI launcher 和 parallel HDF5 必须来自同一套 MPI/编译器环境。例如，用 Open MPI 的 `mpic++` 链接由 MPICH 编译的 HDF5，即使 CMake 能找到文件，也不是可用环境。集群上还应避免同时启用系统库、Conda MPI 和管理员 module 中的 HDF5。

依赖检查应从干净的登录环境开始，并且只加载一套编译器/MPI/HDF5 module。`CMAKE_PREFIX_PATH`、`HDF5_ROOT`、`PKG_CONFIG_PATH` 和 `SDDS_ROOT` 用于明确选择管理员提供的非系统安装。额外 CMake 参数可通过重复指定 `--cmake-arg` 传入安装脚本。

## 无管理员权限的本地依赖引导

联网机器建议直接从远程仓库取得一个全新源码树，避免旧压缩包中的截断文件、CRLF
换行或半成品运行库混入安装：

```bash
git clone --branch <分支或标签> --single-branch <仓库URL> APRL
cd APRL
```

更新已有源码树时使用 `git pull --ff-only`；不要把另一台机器生成的
`aprl-runtime/` 提交或复制进源码仓库。运行库会在目标机器上按目标 MPI/编译器栈
重新生成。

如果机器能够访问互联网，最简单的入口是：

```bash
./install_online.sh --profile core
```

该脚本记录操作系统、架构、调度环境、MPI wrapper 与 module 信息，拒绝在托管集群上
误用 Conda MPI，随后调用下述固定版本、带 SHA-256 校验的依赖引导。它在没有显式选择
时会尝试银杏平台公开的 `intel/openmpi/5.0.3/gcc8.5.0` module；其他集群可使用：

```bash
./install_online.sh --mpi-module <站点提供的MPI模块名> --profile core
```

默认在登录节点跳过两 rank 启动和 MPI 回归测试，但仍执行 APRL 的其余轻量
回归测试；在单节点且至少分配了两个 task 的 Slurm 作业内则默认启用 MPI 检查。单 task
或多节点作业不会擅自使用节点本地临时目录做检查；需要强制检查时可显式传入
`--mpi-run`。正式运行前应在
计算节点分配中进行 MPI 预检。使用 `--plan-only` 可只查看识别结果而不下载或编译。
登录节点默认最多使用 8 个编译进程，allocation 内最多自动使用 16 个；可用
`--jobs N` 明确覆盖。环境摘要保存在 `.cache/aprl-install/`，安装结果和依赖锁位于
`aprl-runtime/`。

自动下载只替代通用编译和运行依赖；多节点 MPI 刻意复用平台 module，因为调度器、
PMI/PMIx 与高速互联支持属于站点环境，不能由普通便携 MPI 安全替代。
例如银杏平台的作业脚本应在运行 APRL 前保留：

```bash
module load intel/openmpi/5.0.3/gcc8.5.0
source /shared/path/APRL/aprl-runtime/share/aprl/runtime-env.sh
```

第一次正式计算前，在一个计算节点 allocation 内执行一次双 rank 预检：

```bash
salloc -N 1 -n 2 -t 00:10:00
module purge
module load intel/openmpi/5.0.3/gcc8.5.0
source "$PWD/aprl-runtime/share/aprl/runtime-env.sh"
./check_dependencies.sh --profile core \
  --mpi-cxx "$(command -v mpic++)" \
  --mpi-launcher "$(command -v mpirun)" \
  --report dependency-report-compute-node.log
```

登录节点的一键安装成功只说明编译、单 rank 链接/运行和非 MPI 回归通过；上述
allocation 内检查才确认两 rank collective parallel-HDF5 与站点启动链路可用。

联网机器可以把经过 SHA-256 校验的固定版本依赖编译到 APRL 安装前缀，不修改操作系统：

```bash
./install.sh --bootstrap-deps --profile core \
  --prefix "$PWD/aprl-runtime" --jobs 8
```

该命令自动准备 CMake（主机缺少合格版本时）、yaml-cpp、zlib 和由当前 MPI C wrapper 编译且启用 deflate 的 parallel HDF5；`full` 配置还准备 FFTW3 double precision 与线程库。源码包缓存于 `.cache/aprl-dependencies`，后续重复安装不重新下载。已经填充缓存的离线机器可增加 `--offline-deps`。依赖标记同时包含源码校验值和当前 MPI/C/C++ 工具链指纹；更换 module 或编译器后不会静默复用旧 ABI 的库。

安装前缀包含 `bin/`、`lib/`、头文件、依赖锁文件和 `share/aprl/runtime-env.sh`。APRL 可执行文件写入相对 `$ORIGIN/../lib`/`lib64` RPATH，因此本地 HDF5、zlib 与 yaml-cpp 不需要安装到每个计算节点；前缀只需位于所有节点可见的共享文件系统。平台 MPI 仍由 module 提供，因此顺序必须是先加载与编译时相同的 module，再载入私有运行环境：

```bash
module load <与编译时相同的MPI模块>
source /shared/path/aprl-runtime/share/aprl/runtime-env.sh
```

托管集群默认复用站点 MPI，并只在本地构建其余库：

```bash
./install.sh --bootstrap-deps --profile core \
  --prefix "$PWD/aprl-runtime" \
  --mpi-cc "$(command -v mpicc)" \
  --mpi-cxx "$(command -v mpic++)" \
  --mpi-launcher mpirun --skip-mpi-run --jobs 8
```

登录节点禁止 MPI 启动时先使用 `--skip-mpi-run`，随后必须在真实调度分配中重新执行 `check_dependencies.sh` 的双 rank 检查。平台 MPI/module 在运行作业中仍须加载，因为它负责 Slurm、PMI/PMIx 和高速互联。普通工作站或不依赖站点网络的单节点环境可显式使用私有 MPICH：

```bash
./install.sh --bootstrap-deps --bundle-mpi --profile core \
  --prefix "$PWD/aprl-runtime" --jobs 8
```

不要在同一前缀内从私有 MPICH 切换为平台 MPI；安装器会拒绝这种混用，应选择新前缀。`glibc`、Linux 内核、基础 C/C++ 编译器、`make`、`tar`、`realpath` 以及 `curl`/`wget` 和 SHA-256 校验工具仍由主机提供；这些是无法安全打包替代的基础宿主能力。私有 MPI 使用通用 TCP socket 通信，是方便的单节点兜底，不等价于针对超算高速网络优化的站点 MPI。

只查看下载、版本和构建计划而不修改文件：

```bash
./bootstrap_dependencies.sh --profile core --plan-only
```

网络必须经过站点镜像时，可成对设置 `APRL_CMAKE_*`、`APRL_YAML_CPP_*`、
`APRL_ZLIB_*`、`APRL_HDF5_*`、`APRL_FFTW_*` 的 `URL` 与 `SHA256` 环境变量。
脚本不会接受未提供匹配 SHA-256 的替代源码，也不会因为使用镜像而关闭校验。

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
Open MPI 5 在 Slurm 分配内优先推荐 `mpirun`。只有站点已明确验证
PMIx direct launch 时才改用 `--mpi-launcher srun`；指定的运行器必须支持通用的
`-n 2` rank 参数。

`--allow-serial-hdf5` 仅用于单 rank 本地调试。默认配置会拒绝 serial HDF5，防止程序成功编译后才在超算多 rank 读粒子时失败。

## 自动安装系统依赖

在可联网、允许包管理器访问的软件机上，以管理员身份执行：

```bash
sudo ./install_dependencies.sh --profile full \
  --report dependency-install.log
```

工具支持 `apt-get`、`dnf`、`yum`、`zypper` 和 `pacman`，安装编译器、CMake、MPI、parallel HDF5、yaml-cpp、FFTW 及 SDDS 所需的系统支持库。包名随发行版和软件源可能变化；若安装失败，完整的包管理器输出仍保留在日志中。

脚本要求以 `root` 身份运行，不会在内部自行调用 `sudo`，也不会修改源码。脚本不安装官方 SDDS。安装结束后以 root-safe 模式重新编译和链接依赖探针，但跳过两 rank 启动；MPI 运行检查仍须在实际作业环境中以普通用户身份完成。

离线集群环境不适用自动包安装器。此类环境应运行 `check_dependencies.sh` 生成需求日志，并由管理员提供一致的 module。具有离线软件源的站点可由管理员使用 `--skip-refresh` 调用安装器。

## 安装程序

默认安装完整配置，包括主程序、通用小工具和五个后处理器：

```bash
./install.sh --prefix /path/to/install-prefix
```

若需要系统安装：

```bash
sudo ./install.sh --prefix /opt/aprl
```

如果需要原生 SDDS 转换器：

```bash
./install.sh --prefix /path/to/install-prefix \
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
- `--mpi-launcher /path/to/mpirun`：明确两 rank 依赖检查使用的运行器；
- `--generator NAME`：明确 CMake generator；
- `--cmake-arg ARG`：传递站点 toolchain 或依赖路径，可重复；
- `--allow-serial-hdf5`：只构建单 rank 本地版本，并跳过 parallel-HDF5 专项测试。

默认前缀是 `/usr/local`，默认一次性构建目录是源码下的 `build-install`。安装完成后，二进制位于 `<prefix>/bin`，文档、示例和绝对路径安装回执位于 `<prefix>/share/aprl`。其中主程序为 `<prefix>/bin/aprl`。

重复安装到同一前缀时，脚本会先保存并验证旧回执；新安装及新回执成功后，才删除旧版本中已经不存在的文件。这样更新不会因某个已移除的工具留下不可追踪的旧二进制。

## 卸载

先检查目标：

```bash
./uninstall.sh --prefix /path/to/install-prefix --dry-run
```

确认后卸载：

```bash
./uninstall.sh --prefix /path/to/install-prefix
```

卸载器只读取 `install.sh` 创建的绝对路径回执，并拒绝任何不在指定前缀内的条目。它删除安装副本和默认的 `build-install` 编译目录，但不删除输出数据、依赖库、系统包或源码。没有合法回执时不会根据文件名猜测并删除内容。

系统前缀下卸载需要相同权限：

```bash
sudo ./uninstall.sh --prefix /opt/aprl
```

自动依赖安装器安装的是系统共享软件包，因此普通程序卸载不会反向删除它们；这些库可能仍被其他程序使用。
