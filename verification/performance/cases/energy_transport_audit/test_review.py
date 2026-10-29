#!/usr/bin/env python3
"""Adversarial accounting-reader checks; no PIC simulation or large input."""
from pathlib import Path
import tempfile
import unittest
import h5py
import numpy as np
from review import energy_review, plane_review


class EnergyReader(unittest.TestCase):
    def fixture(self, path):
        names = ('time_box_s', 'particle_kinetic_active_J', 'particle_kinetic_removed_J', 'field_energy_interior_J',
            'outward_field_energy_J', 'prescribed_work_J', 'boundary_source_work_J', 'balance_residual_raw_J',
            'balance_residual_J', 'moving_inserted_field_J', 'moving_cropped_field_J', 'moving_remap_field_change_J',
            'moving_remap_interface_correction_J', 'inlet_overwrite_work_J', 'longitudinal_boundary_outward_J',
            'frozen_constraint_work_J', 'carrier_constraint_work_J', 'relative_balance_exchange')
        rows = np.zeros(3, dtype=[('step','u8'),('closure_valid','u1'),('removed_macroparticles','u8')]+
            [(name,'f8') for name in names]+[('outward_field_energy_face_J','f8',(6,))])
        rows['step']=[0,1,2]; rows['removed_macroparticles']=[0,0,1]
        rows['time_box_s']=[0,.1,.2]
        rows['particle_kinetic_active_J']=[10,9.9,9.8]
        rows['particle_kinetic_removed_J']=[0,.1,.2]
        rows['field_energy_interior_J']=[1,1.15,1.1]
        rows['prescribed_work_J']=[0,.1,.2]
        rows['outward_field_energy_face_J'][:,0]=[0,.01,.02]
        rows['outward_field_energy_face_J'][:,4]=[0,-.02,-.03]
        rows['outward_field_energy_face_J'][:,5]=[0,.03,.04]
        rows['outward_field_energy_J']=rows['outward_field_energy_face_J'].sum(axis=1)
        rows['moving_inserted_field_J']=[0,.02,.04]
        rows['moving_cropped_field_J']=[0,.01,.02]
        rows['moving_remap_interface_correction_J']=[0,.003,.004]
        rows['moving_remap_field_change_J']=rows['moving_inserted_field_J']-rows['moving_cropped_field_J']+rows['moving_remap_interface_correction_J']
        rows['inlet_overwrite_work_J']=[0,.01,.02]
        rows['frozen_constraint_work_J']=[0,.005,.007]
        rows['carrier_constraint_work_J']=[0,.002,.004]
        rows['longitudinal_boundary_outward_J']=rows['outward_field_energy_face_J'][:,4:].sum(axis=1)-rows['moving_remap_field_change_J']-rows['inlet_overwrite_work_J']
        rows['balance_residual_raw_J']=rows['field_energy_interior_J']-1+rows['outward_field_energy_J']-rows['prescribed_work_J']
        rows['balance_residual_J']=rows['balance_residual_raw_J']-rows['moving_remap_field_change_J']-rows['inlet_overwrite_work_J']-rows['frozen_constraint_work_J']-rows['carrier_constraint_work_J']
        denominator=abs(rows['field_energy_interior_J']-1)+abs(rows['outward_field_energy_J'])+abs(rows['prescribed_work_J'])+abs(rows['moving_remap_field_change_J'])+abs(rows['inlet_overwrite_work_J'])+abs(rows['frozen_constraint_work_J'])+abs(rows['carrier_constraint_work_J'])
        rows['relative_balance_exchange']=np.divide(rows['balance_residual_J'],denominator,out=np.zeros(3),where=denominator>0)
        with h5py.File(path,'w') as stream:
            group=stream.create_group('energy_ledger')
            group.attrs['frame']='boosted computational frame'
            group.attrs['format_version']=4
            group.create_dataset('records',data=rows)

    def test_consistent_nonzero_residual_is_reported_not_hidden(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.fixture(path)
            result=energy_review(path,False,export_csv=False)
            self.assertTrue(result['checks_pass'])
            self.assertEqual(result['closure_valid_values'],[0])
            self.assertNotEqual(result['residual_decomposition']['corrected_final_J'],0)
            self.assertIn('NOT_ASSIGNED',result['conservation_acceptance'])
            self.assertFalse(path.with_suffix('.csv').exists())

    def test_missing_overwrite_accounting_fails_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.fixture(path)
            with h5py.File(path,'r+') as stream:
                rows=stream['energy_ledger/records'][:]
                rows['balance_residual_J']+=rows['inlet_overwrite_work_J']
                stream['energy_ledger/records'][:]=rows
            result=energy_review(path,False,export_csv=False)
            self.assertFalse(result['checks_pass'])
            self.assertFalse(result['identities']['corrected_residual_reconstructed']['pass'])

    def test_wrong_remap_sign_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.fixture(path)
            with h5py.File(path,'r+') as stream:
                rows=stream['energy_ledger/records'][:]
                rows['moving_remap_interface_correction_J']*=-1
                stream['energy_ledger/records'][:]=rows
            result=energy_review(path,False,export_csv=False)
            self.assertFalse(result['checks_pass'])
            self.assertFalse(result['identities']['remap_equals_insert_minus_crop_plus_interface']['pass'])




    def candidate(self, path):
        self.fixture(path)
        fields = ('particle_kinetic_box_removed_event_J','particle_removal_interpolation_delta_J','particle_grid_work_J','prescribed_boris_work_J','magnetic_rotation_delta_K_J',
            'boris_push_delta_K_J','boris_roundoff_residual_J','particle_work_balance_residual_J',
            'field_exchange_residual_J','balance_residual_boris_J','particle_kinetic_lab_active_J',
            'particle_kinetic_lab_removed_J','particle_cPz_box_active_J','particle_cPz_box_removed_J',
            'delta_particle_kinetic_box_J','delta_particle_kinetic_lab_J','delta_particle_cPz_box_J',
            'particle_lorentz_delta_residual_J','delta_field_energy_J','boost_gamma','boost_beta')
        with h5py.File(path,'r+') as f:
            g=f['energy_ledger']; old=g['records'][:]
            rows=np.zeros(len(old), dtype=old.dtype.descr+[(n,'f8') for n in fields])
            for n in old.dtype.names: rows[n]=old[n]
            rows['boost_gamma']=2; rows['boost_beta']=np.sqrt(.75)
            rows['particle_kinetic_box_removed_event_J']=rows['particle_kinetic_removed_J']
            rows['prescribed_boris_work_J']=2*rows['prescribed_work_J']
            rows['particle_grid_work_J']=-rows['prescribed_boris_work_J']
            rows['delta_field_energy_J']=rows['field_energy_interior_J']-rows['field_energy_interior_J'][0]
            corrected=rows['balance_residual_J']+rows['prescribed_work_J']-rows['prescribed_boris_work_J']
            rows['field_exchange_residual_J']=corrected;rows['balance_residual_boris_J']=corrected
            del g['records'];g.create_dataset('records',data=rows);g.attrs['format_version']=5

    def test_candidate_work_split_and_nonclosure_report(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.candidate(path)
            result=energy_review(path,True,export_csv=False)
            self.assertTrue(result['checks_pass'])
            self.assertNotEqual(result['physical_imbalance']['boris']['final_J'],0)
            self.assertIsNone(result['physical_imbalance']['boris']['final_over_abs_net_particle_change'])

    def test_wrong_boris_work_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.candidate(path)
            with h5py.File(path,'r+') as f:
                rows=f['energy_ledger/records'][:];rows['prescribed_boris_work_J']*=2
                f['energy_ledger/records'][:]=rows
            self.assertFalse(energy_review(path,True,export_csv=False)['checks_pass'])


class PlaneReader(unittest.TestCase):
    def fixture(self,path):
        fields=('time_box_s','time_lab_s','first_field_time_lab_s','last_field_time_lab_s',
                'field_power_lab_W','physical_field_power_lab_W','field_energy_lab_J','physical_field_energy_lab_J',
                'field_energy_per_length_lab_J_per_m','physical_field_energy_per_length_lab_J_per_m')
        dtype=[(n,'f8') for n in fields]+[('field_available','u1')]+[(n,'u8') for n in
               ('field_sample_count','outside_field_window_steps','field_gap_count')]
        for prefix in ('particle_','physical_particle_'):
            dtype += [(prefix+n,'u8') for n in ('downstream_count','upstream_count')]+[(prefix+'charge_C','f8')]
            for n in ('kinetic_lab_J','total_lab_J','cpz_lab_J','kinetic_box_J','total_box_J','cpz_box_J'):
                dtype += [(prefix+n,'f8'),(prefix+n+'_low','f8')]
        rows=np.zeros(3,dtype=dtype); count=np.arange(3); rest=count*9.1093837139e-31*299792458**2
        rows['time_box_s']=rows['time_lab_s']=rows['last_field_time_lab_s']=[0,.5,1]
        rows['field_sample_count']=[1,2,3]; rows['field_available']=1
        for prefix in ('','physical_'):
            rows[prefix+'field_power_lab_W']=2;rows[prefix+'field_energy_lab_J']=[0,1,2]
            rows[prefix+'field_energy_per_length_lab_J_per_m']=1
        for prefix in ('particle_','physical_particle_'):
            rows[prefix+'downstream_count']=count;rows[prefix+'charge_C']=-count*1.602176634e-19
            rows[prefix+'total_box_J']=rest;rows[prefix+'total_lab_J']=2*rest
            rows[prefix+'kinetic_lab_J']=rest;rows[prefix+'cpz_lab_J']=2*np.sqrt(.75)*rest
        with h5py.File(path,'w') as f:
            g=f.create_group('energy_plane');g.create_dataset('records',data=rows)
            g.create_dataset('complete',data=1);g.create_dataset('committed_records',data=3)
            g.attrs.update(format_version=1,boost_gamma=2.,boost_beta=np.sqrt(.75))
        return {'particles':2},{'mesh':{'boost_gamma':2.},'beam':{'input':{'electrons':2.}}}

    def test_compensated_particle_four_momentum_and_flux(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'energy.h5';plan,card=self.fixture(path)
            result=plane_review(path,plan,card)
            self.assertTrue(result['checks_pass'])
            self.assertFalse(result['coverage']['complete_means_full_field_history'])

    def test_bad_flux_integral_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'energy.h5';plan,card=self.fixture(path)
            with h5py.File(path,'r+') as f:
                rows=f['energy_plane/records'][:];rows['field_energy_lab_J'][-1]+=1
                f['energy_plane/records'][:]=rows
            self.assertFalse(plane_review(path,plan,card)['checks_pass'])

    def test_complete_with_missing_particles_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'energy.h5';plan,card=self.fixture(path)
            plan['particles']=3
            with self.assertRaisesRegex(RuntimeError,'cohort'): plane_review(path,plan,card)

    def test_unavailable_nan_is_not_a_finite_field(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'energy.h5';plan,card=self.fixture(path)
            with h5py.File(path,'r+') as f:
                rows=f['energy_plane/records'][:];rows['field_power_lab_W'][1]=np.nan
                f['energy_plane/records'][:]=rows
            with self.assertRaisesRegex(RuntimeError,'Nonfinite'): plane_review(path,plan,card)


if __name__=='__main__':
    unittest.main()
