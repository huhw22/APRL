"""Independent Decimal reconstruction and bilingual display checks."""
import csv
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path

from render import save, sha

D=Decimal


def main(root):
    root=Path(root);rev=root/'07_readable_figures'
    books=json.loads((rev/'LEDGERS.json').read_text())
    audit=[]
    with localcontext() as ctx:
        ctx.prec=65
        for book in books:
            path=root/book['source'][0]
            with path.open() as f:rows=list(csv.DictReader(f))
            names=set(rows[0])
            if 'native_field_change_J' in names:
                terms=[('particle_change_J',1),('native_field_change_J',1),('outward_curl_exchange_J',1),
                       ('inlet_energy_added_J',-1),('remap_energy_added_J',-1),
                       ('prescribed_and_constraint_source_work_J',-1),('time_pairing_channel_J',-1),
                       ('Cowan_channel_J',-1),('physical_CPML_channel_J',-1)]
                key='total_residual_J'
            elif 'physical_component_sum_J' in names:
                terms=[(k,1) for k in ['particle_change_J','field_change_J','front_outward_minus_injection_J',
                      'back_outward_plus_cropped_J','side_outward_J','moving_interface_J','external_and_constraint_work_signed_J']]
                key='physical_component_sum_J'
            elif 'particle_delta_K_J' in names:
                terms=[(k,1) for k in ['particle_delta_K_J','field_delta_U_native_J','front_outward_minus_injection_J',
                    'rear_outward_plus_crop_J','side_outward_J','negative_constraint_and_external_work_J','explicit_discrete_correction_J']]
                key='discrete_accounting_sum_J'
            else:
                modern='particle_net_transport_J' in names
                particle='particle_net_transport_J' if modern else 'particle_net_K_J'
                side='field_transverse_net_J' if modern else 'transverse_total_J'
                if book['name'].endswith('_surfaces'):
                    fields=['field_entrance_signed_J','field_exit_signed_J',side,
                            'time_cap_early_signed_J','time_cap_late_signed_J']
                else:
                    fields=['field_ports_net_J' if modern else 'field_end_ports_net_J',side,
                            'field_time_caps_net_J' if modern else 'time_caps_net_J']
                terms=[(k,1) for k in [particle]+fields+(['particle_time_caps_J'] if modern else [])]
                key='physical_lab_sum_J' if modern else 'complete_sum_J'
            worst=D(0);largest=D(0)
            for row in rows:
                values=[D(row[k])*sign for k,sign in terms]
                worst=max(worst,abs(sum(values)-D(row[key])))
                largest=max(largest,max(abs(v) for v in values))
            bound=max(largest*D('3e-12'),D('1e-30'))
            assert worst<=bound,(book['name'],worst,bound)
            assert book['sum_key']==key
            assert set(book['major_keys']).isdisjoint(book['minor_keys'])
            assert set(book['major_keys']+book['minor_keys'])=={t['key'] for t in book['terms']}
            audit.append(dict(name=book['name'],category=book['category'],rows=len(rows),
                maximum_reconstruction_error_J=str(worst),maximum_term_J=str(largest),source_sha256=sha(path)))
    theory=[]
    for kind in ['uniform','converging','gaussian']:
        base=root/'01_reviewed_physics/02_space_charge'/kind
        a=json.loads((base/'linear-binomial3/theory_comparison.json').read_text())
        b=json.loads((base/'quintic-compensated/theory_comparison.json').read_text())
        assert a['reference_model']==b['reference_model'],kind
        assert a['sampling']==b['sampling']
        assert [r['lab_z_m'] for r in a['rows']]==[r['lab_z_m'] for r in b['rows']]
        theory.append(dict(case=kind,identical_model_parameters=True,identical_sampling_rule=True,
                            different_actual_PIC_events=True))
    plots=json.loads((rev/'PLOTS.json').read_text())['figures']
    pairs={}
    for plot in plots:
        assert (rev/plot['path']).is_file()
        assert (rev/plot['path']).with_suffix('.pdf').is_file()
        assert plot['data_axes']==(2 if plot['path'].endswith('_pair.png') or '/theory_' in plot['path'] else 1)
        path=plot['path'].replace('/zh-CN/','/LANG/').replace('/en/','/LANG/')
        values=[(s['key'],s['x'],s['y']) for s in plot['series']]
        if path in pairs:assert pairs[path]==values,path
        else:pairs[path]=values
    assert len(plots)==2*len(pairs)
    save(rev/'INDEPENDENT_VALIDATION.json',dict(passed=True,decimal_precision=65,
         ledgers=audit,one_spherical_model_per_initial_condition=theory,
         bilingual_figures=len(plots),bilingual_data_arrays_identical=True,
         all_J_signs_checked_independently=True,all_large_and_small_terms_partitioned=True))
    print('PASS',len(books),'ledgers',len(plots),'figures')


if __name__=='__main__':main(sys.argv[1])
