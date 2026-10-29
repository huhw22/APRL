"""Normalization edge cases: zero crossings, no particle work, and scale invariance."""
import unittest
import numpy as np
from relative_ledger import construct

def fixture(k,field):
    k=np.asarray(k,dtype=np.longdouble);field=np.asarray(field,dtype=np.longdouble);n=len(k)
    names=['step','delta_particle_kinetic_box_J','prescribed_boris_work_J','boundary_source_work_J','frozen_constraint_work_J','carrier_constraint_work_J']
    rows=np.zeros(n,dtype=[(x,'i8' if x=='step' else np.longdouble) for x in names]);rows['step']=np.arange(n);rows['delta_particle_kinetic_box_J']=k
    cols=['discrete_curl_exchange_J','time_electric_J','time_magnetic_J','current_time_pairing_J','cowan_operator_J','magnetic_cpml_J','electric_cpml_J','inlet_overwrite_J','remap_J']
    sums={x:np.zeros(n-1,dtype=np.longdouble) for x in cols};sums['native_begin_J']=field[:-1];sums['native_end_J']=field[1:]
    r=k+field-field[0]
    q={'discrete_total_residual_J':r,'field_operator_residual_J':r,'particle_grid_coupling_gap_J':np.zeros(n,dtype=np.longdouble)}
    return rows,sums,np.arange(n),q

class RelativeLedgerTests(unittest.TestCase):
    def test_crossing_and_zero_final_change(self):
        values,pct,meta=construct(*fixture([0,2,-1,0],[10,8,11,10]))
        self.assertEqual(meta['fixed_denominator_J'],2)
        np.testing.assert_array_equal(pct['particle_change_pct'],[0,100,-50,0])
        self.assertIsNone(meta['final_total_residual_pct_absolute_final_particle_change'])
        np.testing.assert_array_equal(pct['total_residual_pct'],0)
    def test_zero_particle_work_uses_initial_field(self):
        _,_,m=construct(*fixture([0,0,0],[10,10,10]))
        self.assertEqual(m['fixed_denominator_J'],10)
        self.assertEqual(m['denominator_policy'],'initial_native_field_energy_zero_particle_change')
    def test_zero_initial_energy_uses_peak_field(self):
        _,pct,m=construct(*fixture([0,0,0],[0,1,2]))
        self.assertEqual(m['fixed_denominator_J'],2)
        self.assertEqual(m['denominator_policy'],'peak_native_field_energy_zero_initial_field_and_particle_change')
        np.testing.assert_array_equal(pct['total_residual_pct'],[0,50,100])
    def test_all_zero_is_undefined(self):
        with self.assertRaisesRegex(ValueError,'undefined'):construct(*fixture([0,0],[0,0]))
    def test_rescaling_preserves_percentages_and_signs(self):
        _,a,m=construct(*fixture([0,2,-1],[10,8,12]))
        _,b,n=construct(*fixture([0,200,-100],[1000,800,1200]))
        self.assertEqual(n['fixed_denominator_J'],100*m['fixed_denominator_J'])
        for key in a:np.testing.assert_allclose(a[key],b[key],rtol=1e-17,atol=0)

if __name__=='__main__':unittest.main()
