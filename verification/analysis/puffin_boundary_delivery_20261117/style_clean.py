"""Presentation-only title-and-plot style. Never edits the underlying fields."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection,QuadMesh
from common import sha

def clean_figure(fig):
    # Remove colorbar axes, leaving the actual plotted data axes.
    for ax in list(fig.axes):
        if getattr(ax,'_colorbar',None) is not None:fig.delaxes(ax)
    for ax in fig.axes:
        image_like=bool(ax.images) or any(isinstance(c,QuadMesh) for c in ax.collections)
        if image_like:
            for line in list(ax.lines):line.remove()
        for obj in list(ax.texts):obj.remove()
        for obj in list(ax.patches):obj.remove()
        for obj in list(ax.collections):
            if isinstance(obj,PathCollection):obj.remove()
        if ax.get_legend() is not None:ax.get_legend().remove()
        ax.grid(False);ax.set_axis_off();ax.set_title(ax.get_title(),pad=16)
    title=getattr(fig,'_suptitle',None)
    for txt in list(fig.texts):
        if txt is not title:txt.remove()
    fig.set_size_inches(12.8,7.2,forward=True)
    if len(fig.axes)==2:
        fig.subplots_adjust(left=.035,right=.965,bottom=.035,top=.84,wspace=.08)
        if title is not None:title.set_y(.965)
    else:fig.subplots_adjust(left=.04,right=.96,bottom=.04,top=.88)
    for ax in fig.axes:
        assert not ax.axison and not ax.texts and not ax.patches
        assert ax.get_legend() is None
    assert all(txt is title for txt in fig.texts)

def render_background(root):
    root=Path(root);bg=root/'background';src=bg/'data'
    model=json.loads((bg/'SOURCE_OPERATOR_RESULT.json').read_text())
    m=np.genfromtxt(src/'native_transverse_mask.csv',names=True,delimiter=',')['mask_native'].reshape(81,81)
    tail=np.genfromtxt(src/'tail_mask_comparison.csv',names=True,delimiter=',')
    x=(np.arange(81)-40)*40
    scale=json.loads((bg/'evidence/scale.json').read_text());g=model['grid']
    k=np.linspace(0,60000,301)
    d=np.exp(-g['h_bar']*k*scale['transverse_scale_m']/np.sqrt(2*g['carrier_kz2_bar']))
    records=[]
    for lang in ['zh-CN','en']:
        zh=lang=='zh-CN';out=bg/'figures'/lang;out.mkdir(parents=True,exist_ok=True)
        titles=['横向边界掩膜','横向掩膜剖面','边缘分支谱系数','纵向掩膜几何参考','纵向掩膜对网格长度的敏感性'] if zh else [
          'Transverse boundary mask','Transverse mask profile','Boundary-branch spectral factors','Longitudinal mask and geometric reference','Longitudinal mask sensitivity to grid length']
        names=['01_transverse_mask','02_mask_profile','03_modal_factors','04_longitudinal_length_issue','05_longitudinal_grid_sensitivity']
        for i,name in enumerate(names):
            fig,ax=plt.subplots(figsize=(12.8,7.2));ax.set_title(titles[i])
            if i==0:ax.imshow(m,origin='lower',extent=[-1600,1600,-1600,1600],vmin=0,vmax=1,cmap='viridis',aspect='equal',interpolation='nearest')
            elif i==1:ax.plot(x,m[40],color='#1769a0',lw=2);ax.set_ylim(-.05,1.1)
            elif i==2:
                ax.plot(k,d,color='#1769a0',lw=2);ax.plot(k,d*d,color='#cf661f',lw=2,ls='--');ax.set_ylim(.7,1.025)
            elif i==3:
                a=tail[tail['nz']==96000]
                ax.plot(a['distance_from_tail_um'],a['mask_native'],color='#cf661f',lw=2)
                ax.plot(a['distance_from_tail_um'],a['geometric_parentheses_reference'],color='#1769a0',lw=2,ls='--');ax.set_ylim(-.04,1.08)
            else:
                for n,col,sty in [(95980,'#308367',':'),(96000,'#cf661f','-'),(96020,'#1769a0','--')]:
                    a=tail[tail['nz']==n];ax.plot(a['distance_from_tail_um'],a['mask_native'],color=col,lw=2,ls=sty)
                ax.set_ylim(-.04,1.08)
            clean_figure(fig)
            fig.savefig(out/(name+'.png'),dpi=150);fig.savefig(out/(name+'.pdf'));plt.close(fig)
            records.append(dict(path=str((out/(name+'.png')).relative_to(root)),sha256=sha(out/(name+'.png'))))
    return records
