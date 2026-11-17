"""Copy the selected delivery's unchanged data for a title-only visual revision."""
import json
import shutil
import sys
from pathlib import Path
from common import SOURCE,sha,save

PREVIOUS=SOURCE.parent/'puffin_boundary_explanation_20261117_r3'

def main(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=False)
    old=json.loads((PREVIOUS/'PACKAGE.json').read_text())
    assert sha(PREVIOUS.with_suffix('.zip'))==old['sha256']
    excluded={'PACKAGE.json','SHA256SUMS','RENDER_VALIDATION.json','VALIDATION.json','VISUAL_REVIEW.json','SOURCE_COMMIT.json','figures'}
    for p in PREVIOUS.iterdir():
        if p.name in excluded:continue
        if p.is_dir():
            shutil.copytree(p,root/p.name,ignore=shutil.ignore_patterns('visual_review','figures'))
        else:shutil.copyfile(p,root/p.name)
    selection=json.loads((root/'SELECTION.json').read_text())
    selection['source_delivery']=str(PREVIOUS);selection['source_zip_sha256']=old['sha256']
    save(root/'SELECTION.json',selection)
    save(root/'PRESENTATION.json',dict(title_and_plot_only=True,
       axes_ticks_and_labels=False,colorbars=False,legends=False,
       annotations_and_arrows=False,region_overlays=False,footer_and_statistics=False,
       fixed_intensity_ranges=dict(linear=[0,1],log=[1e-8,1]),
       original_scientific_data_unchanged=True))
    for name in ['README.md','PPT_HANDOFF.md','SUMMARY_EN.md']:
        p=root/name;text=p.read_text()
        text=text.replace('中央方框为 x/y=±520 µm 观察区，原生小域中的竖虚线标记 x方向边缘处理区。每帧增加的是真实传播距离，不是人为移动输入位置。',
          '中央观察区在统计中定义为 x/y=±520 µm。净图版不绘制区域框、坐标、色标或小字；每帧增加的是真实传播距离。强度映射保持固定：线性0–1，对数1e-8–1，均相对同一入射峰值。')
        text=text.replace('两张图在 z=1.505 m 使用同一色标和相同中央方框。','两张图在 z=1.505 m 使用同一物理视野和相同固定强度映射。')
        text=text.replace('同一z、同一色标和相同观察方框','同一z、同一视野和相同固定强度映射')
        text=text.replace('双栏1280×760','双栏1280×720')
        text=text.replace('保留所有色标和范围说明','在讲稿中说明固定强度范围，不在画面上添加色标')
        if name=='README.md':text='> 净图版：图与视频只保留主要图像和标题。坐标、单位、曲线含义和数值解释均放在本文，不嵌入画面。\n\n'+text
        elif name=='PPT_HANDOFF.md':
            text='> 所有视觉素材只保留主图和标题。不要在画面上加回说明小字、箭头、区域框、图例、坐标或数值色标；解释放讲稿。曲线中大域参考为绿实线、原生吸收为蓝虚线。\n\n'+text
        p.write_text(text)
    (root/'data/README.md').write_text((root/'data/README.md').read_text()+'\n净图版只调整显示布局；CSV、JSON科学结果和NPZ数组与r3逐文件字节一致。图像中不显示坐标与色标，单位和固定范围见主文及PRESENTATION.json。\n')
    (root/'README.md').write_text((root/'README.md').read_text()+'''\n## 净图的坐标、颜色与单位\n\n视频及完整横向图的横、纵方向分别为x、y，画面边缘为±1620µm，原节点中心为±1600µm。01、02传播图水平方向为x，竖直方向为z=0–2.1m，从下向上增加。10中央误差图只显示x/y=±520µm。视频和横向强度图均相对同一入射峰值；线性映射0–1、对数映射1e-8–1，颜色由暗到亮表示增强。\n\n03–07曲线图水平方向为传播距离z（m）；03、04、05分别为全域、中央及中央后期放大的积分场强占入射全域值的百分比，06为中央复场差L2占入射全域L2的百分比，07为负kx模式平方范数占入射值的份额（对数纵向）。大域参考是绿实线，原生吸收是蓝虚线；06为单条蓝线。精确坐标和值可查CSV，不以无刻度画面代替定量读数。\n''')
    (root/'background/README.md').write_text((root/'background/README.md').read_text()+'\n附录图同样只保留主图和标题。掩膜图范围0–1；支路系数蓝线为D、橙虚线为D²。纵向比较橙线为原生值、蓝虚线为几何参考；长度敏感性绿点线/橙实线/蓝虚线分别对应95980/96000/96020节点。相关坐标及数值在背景CSV中。\n')
    save(root/'metadata/VISUAL_REVISION_SOURCE.json',dict(previous_zip=str(PREVIOUS.with_suffix('.zip')),
        previous_zip_sha256=old['sha256'],data_hashes={str(p.relative_to(PREVIOUS)):sha(p) for p in (PREVIOUS/'data').iterdir() if p.is_file() and p.name!='README.md'},
        summary_sha256=sha(PREVIOUS/'SUMMARY.json')))
    print(root)

if __name__=='__main__':main(sys.argv[1])
