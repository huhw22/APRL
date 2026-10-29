# 空间电荷理论参照

交付目录：工作区
`results/energy_ledger_test_record_20261030/physical_benchmarks/space_charge_reference_20261030/`。
这是独立于 APRL 的理论参考准备，不是已完成的 PIC 验收。

`reference.py` 用冷球壳隐式解生成均匀球／高斯理论曲线，另算相干向内
速度斜率的均匀球。Gauss–Legendre 积分与空间 E² 积分核对自能，
时间 ODE 独立核对隐式轨迹，监测高斯壳层交叉条件。它不调用 APRL。

```bash
python3 source/aprl/verification/analysis/space_charge_reference/reference.py \
  --output /absolute/new/theory-directory --electrons 100000 --sigma-m 0.0001
```

只有 `reference.json` 明确标为理论输出时，`--refresh` 才允许刷新该目录。
示例参数不是生产卡。保留内部运动非相对论和壳层交叉前的适用条件。

`audit_fixture.py` 只读已有固定 lab 面 H5，计算过面事件的位置—速度
相关性和线性斜率；它不把该输入解释成共同时间的冷束团，也不量测实际
能量变化的成因：

```bash
python3 source/aprl/verification/analysis/space_charge_reference/audit_fixture.py \
  --card /absolute/existing/run.yaml --output /absolute/new/existing_fixture_scope.json
```

`theory_template.md` 用结果目录的相对路径引用图片和数据，复制成该目录的
`README.md` 阅读。正文区分文献模型、自行推导、纯理论预测和未来 PIC 比较。
依赖现有 physics Python 环境的 NumPy、SciPy、Matplotlib、h5py、PyYAML。
