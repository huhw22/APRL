# 十万电子同步辐射算例的侧面 CPML 初态

此例复现计算 20104 的完整输入：十万个电子、真实自场反馈、偶极偏转、原六个场屏、Cowan 移动窗、纵向及 x 向三点核和四侧各 16 格 CPML。默认读入 OSIRIS 思路的四侧匹配初态，在初始磁场半步前装入四面和角区的电磁场、八条卷积记忆，以及经过低频续接的纵向均值。备选的 20056 初态使用相同束团和网格，装入左壁初始场和相应记忆，包括靠内侧八排的单记忆修正。

## 准备和运行

准备已安装的 APRL 2.1、系统 MPI，以及发行时单独提供的 `APRL-2.1.0-beam-cpml-100k.tar.gz` 和同名的 `.sha256` 文件。在数据包所在目录执行：

```bash
sha256sum -c APRL-2.1.0-beam-cpml-100k.tar.gz.sha256
tar -xzf APRL-2.1.0-beam-cpml-100k.tar.gz
python3 beam-cpml-100k/prepare.py --data beam-cpml-100k --verify
python3 beam-cpml-100k/prepare.py --data beam-cpml-100k --output run-20104
cd run-20104
aprl --version
export OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close
mpiexec -n 16 aprl run.yaml
```

原计算使用 16 个 MPI 进程，每进程 8 个线程。请使用集群适配的 MPI 启动器，并为此布局分配相应的资源。脚本只准备输入，输出目录需为新路径；程序启动时的当前目录需是含有生成的 `run.yaml` 的目录。输入卡内的文件路径相对于该目录解析。二进制通过安装版的 PATH 选择，场、屏及账本写入本次运行目录。用于准备运行的 Python 脚本只需要 Python 3 标准库。

若选择左壁初态，准备新的目录：

```bash
python3 beam-cpml-100k/prepare.py --data beam-cpml-100k --method left-20056 --output run-left
```

可从已安装的 `share/aprl/examples/beam_cpml_100k/prepare.py` 调用同一命令，`--data` 指向用户自己的解压位置。准备时校验束团文件、所选方法各电磁场和记忆数组的 SHA-256。此输入保留正常粒子受力与自洽电流沉积，并且在 YAML 中显式选择横向平滑、开放横向参考势、初态方法与数据目录。

## 适用范围

这些交错网格数组针对给定粒子文件、初始位置、束流参数、参考域 416×140×2816、移动窗 416×140×1457、初始纵向偏移 623 格、原 CPML 参数和吸收层厚度制作。采用小端 IEEE 754 双精度。切换束团、参考域、网格间距、步长、层宽或启动时刻时需为新状态求解相应场和记忆并校验。随包发放的两个初态都是固定输入。OSIRIS 式初态可减小部分启动瞬态；已测同一束团的近侧壁仍可产生明显的近场响应。这里的边界选择不表示近场污染已经消除，也不把束团工况扩展为任意几何通用吸收。

## 维护者重建数据包

在源码示例目录运行 `package_inputs.py`，分别提供原 20104 束团 HDF5、四侧修正后的数据目录、20056 左侧数据目录，以及此前原生 Cowan 参考目录（含 `metadata.txt`、`rho.bin`）。打包脚本逐项核对冻结的上游 SHA-256，把所有数据作为普通文件写入独立可搬移目录，输出运行时 `startup.yaml`、`SHA256SUMS`、压缩包及外部 SHA-256。原始文件不会写回。此步骤需要 PyYAML、C99 编译器及 gzip 或 pigz。运行样例本身直接使用已归档数组，无需重新求解参考模态。

```bash
python3 package_inputs.py --beam /path/to/beam.h5 --osiris /path/to/four-side-arrays \
  --left-20056 /path/to/left-arrays --reference /path/to/native-reference --output /path/to/release
```
