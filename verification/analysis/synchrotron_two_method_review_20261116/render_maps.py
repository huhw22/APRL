"""Fixed page margins for standalone field maps and their colour scales."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,colors
r=Path(sys.argv[1]);font=sys.argv[2]
font_manager.fontManager.addfont(font);family=font_manager.FontProperties(fname=font).get_name()
plt.rcParams.update({'font.family':[family,'DejaVu Sans'],'font.size':12,'axes.titlesize':16,'axes.labelsize':13,
 'text.parse_math':False,'axes.unicode_minus':False,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
records=[]
for seed in [20261111,20261112]:
 with np.load(r/'radiation'/f'F3_{seed}_maps.npz') as z:
  x,y=z['x'],z['y'];values={k:z[k] for k in ['LW_full','LW_acc']};values['F3']=z['PIC']
 with np.load(r/'radiation'/f'MM_{seed}_maps.npz') as z:values['MM']=z['PIC']
 vmax=max(values[k][10:-10,10:-10].max() for k in ['LW_full','F3','MM'])*1e6
 for lang in ['zh-CN','en']:
  zh=lang=='zh-CN'
  names={'LW_full':'规定轨迹 LW 总场' if zh else 'Prescribed-orbit full LW',
   'LW_acc':'规定轨迹 LW 加速度项' if zh else 'Prescribed-orbit acceleration LW',
   'F3':'自洽三点核总场' if zh else 'Self-consistent F3 total field',
   'MM':'自洽五次补偿总场' if zh else 'Self-consistent quintic total field'}
  for key,title in names.items():
   fig,ax=plt.subplots(figsize=(9.5,6.6));fig.subplots_adjust(left=.12,right=.83,bottom=.14,top=.89)
   im=ax.imshow(values[key][10:-10,10:-10]*1e6,origin='lower',aspect='equal',cmap='magma',
    extent=[(x[10]-2.5e-6)*1e6,(x[-11]+2.5e-6)*1e6,(y[10]-2.5e-6)*1e6,(y[-11]+2.5e-6)*1e6],
    norm=colors.PowerNorm(gamma=.45,vmin=0,vmax=vmax))
   ax.set(title=title+f' 8 eV {seed}',xlabel='原实验室横坐标（µm）' if zh else 'Original laboratory x (µm)',
    ylabel='横向位置 y（µm）' if zh else 'Transverse y (µm)')
   ca=fig.add_axes([.86,.15,.025,.70]);fig.colorbar(im,cax=ca,
    label='谱场强度（µJ / eV / m2）' if zh else 'Spectral field intensity (µJ / eV / m2)')
   fig.canvas.draw()
   for obj in [ax.title,ax.xaxis.label,ax.yaxis.label,ca.yaxis.label]:
    b=obj.get_window_extent(fig.canvas.get_renderer())
    assert b.x0>=0 and b.y0>=0 and b.x1<=fig.bbox.width and b.y1<=fig.bbox.height,obj.get_text()
   for ext in ['png','pdf']:fig.savefig(r/'figures'/lang/f'{seed}_{key}_8eV.{ext}',dpi=170)
   plt.close(fig);records.append(dict(name=f'{seed}_{key}_8eV',language=lang,data_axes=1,color_scale_max=float(vmax),fixed_margins_checked=True))
p=r/'PLOTS.json';old=json.loads(p.read_text());replaced={(x['name'],x['language']) for x in records}
old['figures']=[x for x in old['figures'] if (x['name'],x['language']) not in replaced]+records
assert len(old['figures'])==80;p.write_text(json.dumps(old,ensure_ascii=False,indent=2)+'\n')
print('MAPS_FIXED',len(records),flush=True)
