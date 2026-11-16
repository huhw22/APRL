"""Resolve small internal exchanges after the independently recorded magnet work."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,ticker
r=Path(sys.argv[1]);font=sys.argv[2];font_manager.fontManager.addfont(font)
family=font_manager.FontProperties(fname=font).get_name()
plt.rcParams.update({'font.family':[family,'DejaVu Sans'],'font.size':12,'axes.titlesize':16,'axes.labelsize':13,
 'text.parse_math':False,'axes.unicode_minus':False,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
records=[]
for case in json.loads((r/'RESULT.json').read_text())['cases']:
 a=np.genfromtxt(r/'cases'/case['tag']/'boost/physical_components.csv',delimiter=',',names=True,dtype=np.longdouble)
 effective=a['particle_change_J']+a['external_prescribed_work_signed_J']
 sources=a['external_and_constraint_work_signed_J']-a['external_prescribed_work_signed_J']
 ys=[effective,a['field_change_J'],a['front_outward_minus_injection_J'],a['back_outward_plus_cropped_J'],
     a['side_outward_J'],a['moving_interface_J'],sources,a['physical_component_sum_J']]
 assert np.max(abs(sum(ys[:-1])-ys[-1]))<1e-23
 for lang in ['zh-CN','en']:
  zh=lang=='zh-CN';fig,ax=plt.subplots(figsize=(10.8,6.6),layout='constrained')
  labels=['粒子变化减外磁功','内场能量变化','前端向外减注入','后端向外加裁出','四侧净传出','移窗接口项','约束及边界源功（负号）','物理分项加和'] if zh else [
   'Particle change minus magnetic work','Field change','Front minus injection','Back plus cropped','Transverse outward',
   'Moving interface','Negative constraint and boundary work','Sum']
  ss=[]
  for i,(name,y) in enumerate(zip(labels,ys)):
   x=np.asarray(a['time_box_s']*1e12,float);v=np.asarray(y*1e15,float)
   ax.plot(x,v,label=name,lw=1.8,color='#111111' if i==7 else None,linestyle=['-','--','-.',':'][i%4])
   ss.append(dict(label=name,x=x.tolist(),y=v.tolist()))
  method=('三点核' if case['method']=='F3' else '五次补偿') if zh else case['method']
  kind=('弯铁' if case['kind']=='bend' else '直线') if zh else case['kind']
  ax.set(title=f'{method} {kind} '+('扣除已记外磁功后的交换' if zh else 'after recorded magnetic-work removal'),
   xlabel='计算系模拟时间（ps）' if zh else 'Computational time (ps)',ylabel='带符号能量（fJ）' if zh else 'Signed energy (fJ)')
  ax.axhline(0,color='#888888',lw=.7);ax.grid(alpha=.18);ax.legend(frameon=False,ncol=2,fontsize=9.5)
  ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p:f'{x:g}'))
  name=case['tag']+'_boost_internal_balance'
  for ext in ['png','pdf']:fig.savefig(r/'figures'/lang/(name+'.'+ext),dpi=170)
  plt.close(fig);records.append(dict(name=name,language=lang,data_axes=1,series=ss,recorded_external_work_removed=True))
p=r/'PLOTS.json';old=json.loads(p.read_text());keys={(x['name'],x['language']) for x in records}
old['figures']=[x for x in old['figures'] if (x['name'],x['language']) not in keys]+records
assert len(old['figures'])==96;p.write_text(json.dumps(old,ensure_ascii=False,indent=2)+'\n')
print('INTERNAL_FIGURES',len(records),flush=True)
