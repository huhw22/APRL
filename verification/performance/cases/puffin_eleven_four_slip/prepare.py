"""Freeze templates and an external-field-only source extension; no simulation."""
import math
import shutil
import sys
from pathlib import Path
import numpy as np
from common import W, C, read, save, sha, need


def main(root):
    root=Path(root).resolve(); root.mkdir(exist_ok=False)
    template=Path(__file__).resolve().parent
    shutil.copytree(template, root/'scripts', ignore=shutil.ignore_patterns('__pycache__'))
    cases=template.parent
    for name, source in [('source_math.py', cases/'physical_electron_source/source_math.py'),
                         ('source_verify.py', cases/'physical_electron_source/verify.py'),
                         ('binding.py', cases/'source_coupling_integration/binding.py'),
                         ('puffin_units.py', cases/'puffin_same_sample_20261113/prepare.py'),
                         ('legacy_validate.py', cases/'puffin_same_sample_20261113/validate.py')]:
        shutil.copy2(source, root/'scripts'/name)
    old=W/'results/sparc_sase_11segment_2slip_20261111_r1'
    prior=read(old/'plan.json'); lattice=read(old/'genesis/lattice.json')
    config=read(W/'results/sparc_random_startup_20261112_r2/PLAN.json')['source_config']
    nominal_slip=11*77*prior['physics']['wavelength_m']; duration=4*nominal_slip/C
    config['profile'].update(full_duration_s=duration, window_start_s=0., window_stop_s=duration)
    sys.path.insert(0, str(cases/'physical_electron_source'))
    from source_math import nominal_electrons
    config['electrons']=round(nominal_electrons(config))
    old_puf=W/'results/puffin_same_sample_20261113_r1'
    shutil.copytree(old_puf/'source', root/'source')
    file=root/'source/puffin/lib/undulator/bfields.f90'; original=file.read_text()
    marker="  if (und%undulator_type == 'aprl-ramp77') then"
    need(original.count(marker)==1, 'Magnet insertion point changed')
    addition="  if (und%undulator_type == 'aprl-sparc11') then\n    call getAPRLEleven(sx, sy, sz, bxj, byj, bzj, frame)\n    return\n  end if\n\n"
    original=original.replace(marker,addition+marker)
    phases=', &\n      '.join(f'{p:.17e}_wp' for p in lattice['segment_phases_rad'])
    external=(template/'external_field.f90.in').read_text().replace('@PHASES@',phases)
    pos=original.lower().rfind('end module')
    file.write_text(original[:pos]+external+'\n'+original[pos:])
    shutil.copy2(template/'field_probe.f90', root/'source/field_probe11.f90')
    cmake=root/'source/CMakeLists.txt'
    cmake.write_text(cmake.read_text()+'''\nadd_executable(field_probe11 field_probe11.f90)
target_include_directories(field_probe11 PRIVATE ${CMAKE_BINARY_DIR}/puffin ${HDF5_INCLUDE_DIRS})
target_link_libraries(field_probe11 puffin_mpi)
''')
    build=W/'build'/('puffin-eleven-four-slip-20261113-'+root.name.rsplit('_',1)[-1])
    outputs=[]; dz=.028/80
    for j in range(11):
        nominal=.01+j*(77*.028+.4)+77*.028
        step=math.ceil((nominal-.01)/dz-1e-8)
        outputs.append(dict(segment=j+1,step=step,nominal_lab_z_m=nominal,actual_lab_z_m=.01+step*dz))
    plan=dict(root=str(root),source_config=config,source_seed=config['seed'],target_weight=64,
        sample_records=round(config['electrons']/64),segments=11,periods_per_segment=77,
        nominal_total_slippage_m=nominal_slip,slippage_multiple=4.,duration_s=duration,
        source_charge_C=config['electrons']*1.602176634e-19,source_length_m=C*duration,
        wavelength_m=prior['physics']['wavelength_m'],K=prior['physics']['K'],period_m=.028,
        peak_current_A=config['profile']['peak_current_A'],lattice=lattice,
        binary=str(build/'puffin/puffin'),build=str(build),field_probe=str(build/'field_probe11'),
        upstream_parent=str(old_puf),external_type='aprl-sparc11',propagation_module_periods=990,
        physical_line_length_m=lattice['physical_line_length_m'],integration_tail_m=990*.028-lattice['physical_line_length_m'],
        steps_per_period=80,transverse_nodes=81,dx_dy_m=4e-5,intervals_per_wavelength=20,
        execution=dict(preparation_nodes=1,preparation_ranks=128,main_nodes=2,main_ranks=256),
        physics_timeout_s=46*3600,analysis_timeout_s=3600,outputs=outputs,
        steady_retarded_window_s=[1.5*nominal_slip/C,2.5*nominal_slip/C],
        aperture_half_width_m=.00056,band_half_fraction=.05,
        only_Puffin_submitted=True,APRL_input_saved=True,noise='Independent sample, no quiet clones or extra noise; w64 macro noise retained.',
        loading='Full smooth parabolic support, independent fixed-N physical source uniformly thinned without replacement.',
        magnet_implementation='One continuous module contains all eleven original magnets, ten .4m gaps and finite quadrupoles; no intermediate correctTrans calls.',
        output_position_convention='First integration step at/after each true segment exit, <=.35mm delay; actual coordinates retained for later PIC planes.')
    save(root/'PLAN.json',plan)
    (root/'inputs').mkdir(); (root/'preparation').mkdir()
    for name in ['import-check','main','smoke-matched']: (root/name).mkdir()
    (root/'smoke-stock').symlink_to(old_puf/'smoke-stock',target_is_directory=True)
    dst=root/'smoke-matched/inputs/3D'; dst.mkdir(parents=True)
    for name in ['clara_test.in','beam_file_h5.in','seed_file.in','clara_electrons_0.h5']:
        shutil.copy2(old_puf/'smoke-matched/inputs/3D'/name,dst/name)
    files=[root/'PLAN.json',*[p for p in (root/'scripts').iterdir() if p.is_file()],
           *[p for p in (root/'source').rglob('*') if p.is_file()],*[p for p in dst.iterdir() if p.is_file()]]
    save(root/'FROZEN_PREPARATION.json',{str(p):sha(p) for p in files})
    print(root)


if __name__=='__main__': main(sys.argv[1])
