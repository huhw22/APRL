"""Create a resource-only single-node relaunch without editing frozen runs."""
import copy
import shutil
import sys
from pathlib import Path


def main(campaign,destination):
    campaign=Path(campaign).resolve();destination=Path(destination).resolve()
    sys.path.insert(0,str(campaign/'scripts'))
    from common import read,save,sha,need,frozen
    source=campaign/'eleven'
    need(read(campaign/'PREPARATION_COMPLETE.json')['complete'],'Campaign preparation incomplete')
    need(read(source/'PREPARATION_COMPLETE.json')['complete'],'Eleven-segment admission incomplete')
    need(not (source/'STARTED.json').exists(),'Use a separate restart procedure for a started production')
    frozen(source/'FROZEN_PHYSICS.json')
    destination.mkdir(exist_ok=False)
    shutil.copytree(campaign/'scripts',destination/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
    template=Path(__file__).resolve().parent
    shutil.copy2(template/'single_node.sbatch',destination/'scripts/single_node.sbatch')
    shutil.copy2(Path(__file__),destination/'scripts/relaunch_single_node.py')
    (destination/'main').mkdir();(destination/'inputs').mkdir()
    for path in (source/'main').iterdir():
        need(path.is_file() and path.name in ['main.in','beam.in','matched.latt','writes.in'],'Unexpected original production file')
        shutil.copy2(path,destination/'main'/path.name)
    # All linked data are immutable inputs. Never open either link for writing.
    (destination/'particles.h5').hardlink_to(source/'particles.h5')
    info=read(source/'INPUT_VALIDATED.json');input_path=Path(info['input'])
    (destination/'inputs'/input_path.name).hardlink_to(input_path)
    for name in ['INPUT_VALIDATED.json','NATIVE_VALIDATED.json']:
        shutil.copy2(source/name,destination/name)
    old=read(source/'PLAN.json');plan=copy.deepcopy(old)
    plan.update(root=str(destination),main_nodes=1,main_ranks=128)
    plan['execution'].update(main_nodes=1,main_ranks=128)
    save(destination/'PLAN.json',plan)
    differences=[key for key in old if old[key]!=plan[key]]
    need(set(differences)=={'root','main_nodes','main_ranks','execution'},'Non-resource plan change')
    original_execution=copy.deepcopy(old['execution'])
    original_execution.update(main_nodes=1,main_ranks=128)
    need(plan['execution']==original_execution,'Non-resource execution change')
    unchanged={}
    for name in ['main.in','beam.in','matched.latt','writes.in']:
        need(sha(source/'main'/name)==sha(destination/'main'/name),'Physical deck changed')
        unchanged[name]=sha(destination/'main'/name)
    need(sha(destination/'particles.h5')==info['native_sha256'],'Native events changed')
    need(sha(input_path)==info['input_sha256'],'Independent source changed')
    gate=read(source/'PREPARATION_COMPLETE.json')
    gate.update(reused_from=str(source),resource_only_relaunch=True,actual_import_check_layout='1 node,128MPIx1',
        physics_and_native_inputs_byte_identical=True,new_production_started=False)
    save(destination/'PREPARATION_COMPLETE.json',gate)
    save(destination/'RESOURCE_ONLY_VALIDATION.json',dict(passed=True,original=str(source),
        changed_plan_keys=differences,old_nodes=2,old_ranks=256,new_nodes=1,new_ranks=128,
        physical_decks_sha256=unchanged,input_sha256=info['input_sha256'],native_sha256=info['native_sha256'],
        physical_source_shared=True,immutable_input_hardlinks=True,binary_unchanged=True,
        original_physics_timeout_s=old['physics_timeout_s'],new_physics_timeout_s=plan['physics_timeout_s'],
        reused_native_gate=str(source/'NATIVE_VALIDATED.json'),
        source_config_equal=plan['source_config']==old['source_config'],
        magnetic_lattice_equal=plan['magnetic_lattice']==old['magnetic_lattice'],
        no_new_generation=True,no_new_physical_preflight=True,no_running_frozen_files_modified=True))
    files=[destination/'PLAN.json',destination/'particles.h5',input_path,
        destination/'inputs'/input_path.name,destination/'INPUT_VALIDATED.json',destination/'NATIVE_VALIDATED.json',
        destination/'PREPARATION_COMPLETE.json',destination/'RESOURCE_ONLY_VALIDATION.json',Path(plan['binary']),
        *[p for p in (destination/'scripts').iterdir() if p.is_file()],
        *[p for p in (destination/'main').iterdir() if p.is_file()]]
    save(destination/'FROZEN_PHYSICS.json',{str(p):sha(p) for p in files})
    frozen(source/'FROZEN_PHYSICS.json');frozen(destination/'FROZEN_PHYSICS.json')
    print(destination,flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
