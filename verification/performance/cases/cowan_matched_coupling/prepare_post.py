"""Repair the reader in a separate output tree; retain the original FAILED job."""
import argparse,json,shutil
from pathlib import Path
from prepare import sha,save,require
from prepare_flow import patch_paired_reader
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--flow-root',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();old=a.flow_root.resolve();r=a.output.resolve();require((old/'PIC_COMPLETE.json').exists(),'No closed science marker')
 for name in ['undulator-on','drift-on']:
  text=(old/name/'aprl.log').read_text();require(all(x in text for x in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']),'Science failed to close '+name)
 r.mkdir(parents=True,exist_ok=False);shutil.copytree(old/'frozen_source',r/'frozen_source',ignore=shutil.ignore_patterns('__pycache__','*.pyc'));patch_paired_reader(r/'frozen_source/paired_post.py')
 for name in ['undulator-on','drift-on']:
  d=r/name;d.mkdir()
  for q in (old/name).iterdir():
   if q.is_dir(): (d/q.name).symlink_to(q,target_is_directory=True)
  for n in ['run.yaml','aprl.log']:(d/n).symlink_to(old/name/n)
 (r/'bin').symlink_to(old/'bin',target_is_directory=True);(r/'PIC_COMPLETE.json').symlink_to(old/'PIC_COMPLETE.json');(r/'CALCULATION_COMPLETE.json').symlink_to(old/'CALCULATION_COMPLETE.json');save(r/'PLAN.json',dict(original_flow=str(old),science_replayed=False,original_batch_exit_code=(old/'EXIT_CODE').read_text().strip(),all_raw_directories_read_only_links=True,reader_change='Add independently observed particle_grid_work minus raw applied-current-mid work and original constraint works; original threshold unchanged',old_failed_reader_retained=True));shutil.copy2(Path(__file__).resolve().parent/'post.sbatch',r/'run.sbatch');text=(r/'run.sbatch').read_text().replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r));(r/'run.sbatch').write_text(text);files=[r/'PLAN.json',r/'run.sbatch',*[x for x in (r/'frozen_source').iterdir() if x.is_file()],(r/'bin/replay').resolve()];save(r/'FROZEN.json',{str(x):sha(x) for x in files});print(r/'run.sbatch')
if __name__=='__main__':main()
