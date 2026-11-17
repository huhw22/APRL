# 2.0 发布构建与验收

正式发行必须记录main提交、标签、构建依赖和测试结果，不能只给源码压缩包。

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ -DBUILD_TESTING=ON \
  -DAPRL_BUILD_POSTPROCESSORS=ON -DAPRL_BUILD_SASE_BUNCH_GENERATOR=ON \
  -DAPRL_ENABLE_MOVING_WINDOW=ON -DAPRL_PORTABLE_INSTALL=ON \
  -DCMAKE_INSTALL_PREFIX="$PWD/aprl-runtime"
cmake --build build-release -j4
ctest --test-dir build-release -L required --output-on-failure --timeout 180
cmake --install build-release
./aprl-runtime/bin/aprl --version
./aprl-runtime/bin/aprl --help
```

依赖可能需要CMAKE_PREFIX_PATH/HDF5_ROOT；MPI与parallel HDF5须属于同一栈。
移动窗构建还需FFTW。完整后处理需要FFTW与线程库。SDDS是可选官方库，
缺失时应在交付记录明确标出，不阻断已支持的HDF5工作流。

必需回归覆盖两个算法的配置默认/正选/拒绝旧方法、实际短运行、运行身份、
Maxwell/CPML、连续性、MPI一致性、探测和停止、输出保护以及后处理链。
低层Boris及滤波数学单元测试可以保留，生产读卡没有这些选择。

安装验收需直接调用安装前缀中的aprl，分别运行三点核和五次补偿小卡，
检查manifest版本、实际kernel和完成状态；以旧Boris/linear/linear-notch卡
确认拒绝。执行ldd确认无缺失依赖，重新定位运行前缀后再次验证。
平台MPI/系统库由目标环境提供；包内私有HDF5/yaml等库及依赖说明随运行包提供。

软件回归通过不更新主README中的物理误差结论。六段FEL的lab差额保留，
不能用测试PASS来宣称其已闭合。当前结果图及原始数据不随发行构建重算。
