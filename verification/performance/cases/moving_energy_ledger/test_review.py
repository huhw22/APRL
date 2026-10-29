#!/usr/bin/env python3
"""Adversarial accounting-reader checks; no PIC simulation or large input."""
from pathlib import Path
import tempfile
import unittest
import h5py
import numpy as np
from review import energy_review


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
            group.create_dataset('records',data=rows)

    def test_consistent_nonzero_residual_is_reported_not_hidden(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.fixture(path)
            result=energy_review(path,True,export_csv=False)
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
            result=energy_review(path,True,export_csv=False)
            self.assertFalse(result['checks_pass'])
            self.assertFalse(result['identities']['corrected_residual_reconstructed']['pass'])

    def test_wrong_remap_sign_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ledger.h5';self.fixture(path)
            with h5py.File(path,'r+') as stream:
                rows=stream['energy_ledger/records'][:]
                rows['moving_remap_interface_correction_J']*=-1
                stream['energy_ledger/records'][:]=rows
            result=energy_review(path,True,export_csv=False)
            self.assertFalse(result['checks_pass'])
            self.assertFalse(result['identities']['remap_equals_insert_minus_crop_plus_interface']['pass'])


if __name__=='__main__':
    unittest.main()
