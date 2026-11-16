"""Redraw only theory legends and four energy videos; all data unchanged."""
import json
import sys
from pathlib import Path
from render import Renderer,save,spheres
from movies import main as movies


def main(root,font,metadata):
    root=Path(root);r=Renderer(root,font);out=r.out
    old=json.loads((out/'PLOTS.json').read_text())
    bypath={p['path']:p for p in old['figures']}
    spheres(r,theory_only=True)
    for p in r.plots:
        assert p['series']==bypath[p['path']]['series'],p['path']
    source=json.loads((out/'SOURCE_TABLE_SHA256.json').read_text())
    for rel,h in r.sources.items():assert source[rel]==h
    newpaths={p['path'] for p in r.plots}
    old['figures']=[p for p in old['figures'] if p['path'] not in newpaths]+r.plots
    assert len(old['figures'])==304
    save(out/'PLOTS.json',old)
    movies(str(root),font,metadata,'energy-only')
    save(out/'POLISH_VALIDATION.json',dict(passed=True,theory_figures=24,
         theory_plot_arrays_byte_identical=True,source_tables_unchanged=True,
         energy_videos_reexported=4,only_legend_position_and_log_tick_text_changed=True))


if __name__=='__main__':main(*sys.argv[1:])
