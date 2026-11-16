# 展示修订：清楚地区分平面、理论取样与能量分量

交付内先看 [图片与视频入口](展示入口.md)，再看 [账本和图像逐项说明](能量账本与图像阅读说明.md)。静态图提供 PNG/PDF，视频和图片均保留中文和英文。

此目录从已闭合交付的 CSV 与时间帧 NPZ 重绘，不读取科学 H5、不推进粒子或场、不拟合或强制闭合。原表格、公式及旧版交付保持原值；PPT 素材表指向新版展示图。

账本的两张图合起来包含全部加和分量，另附同一完整加和的独立图。纵向通过、横向侧流和时间端面分别保留。球形对照使用一套相对论球壳模型，按各方法实际过面时刻取样；尺寸下方增加 nm 差值，区分主图几乎重合的曲线。

自由传播视频明确标出 xz/xy 平面及 Ey 分量。四个视图也分别导出；xy 沿 y 的弱场由沿 y 的电 Hertz 矢量角分布解释，详细公式见说明文档。

`LEDGERS.json` 记录每份账本的来源、分项、分组和重建误差；`PLOTS.json` 记录真实绘制数组的校验信息；`THEORY_SAMPLING.json` 记录两种事件取样的差别。视频帧数与时长核验在 `00_vacuum/MOVIE_VALIDATION.json`，总验证在 `VALIDATION.json`。

## English reading key

Use the English figures beside the Chinese figures. A paired budget places major channels on the left and smaller channels plus the complete signed sum on the right. Read each axis unit: both axes use absolute energies and linear scales, with no per-curve normalization. All channels together form one budget. The separate sum figure retains the measured nonzero residual.

The spherical reference is one relativistic shell model, evaluated at each method's actual per-particle crossing events. The black reference line uses F3 events; the narrow grey envelope covers the two event-sampling values. The lower panel shows each simulation minus its own event-matched reference. Energy-change standard deviation is distinct from projected energy spread.

The vacuum color scale shows signed Ey. The electric Hertz vector points along y, suppressing the leading field along that axis in the xy plane. Spherical wavefronts do not imply isotropic vector-field amplitude. The xz and xy cuts share one fixed color scale across all saved times. The analytical movie is an unbounded-vacuum reference and has no absorbing layer.

Computational discrete budgets, physical laboratory budgets, and spectrally filtered radiation diagnostics have different definitions and domains. Their residuals are not interchangeable. All original scientific CSV/JSON/NPZ records remain byte-identical.
