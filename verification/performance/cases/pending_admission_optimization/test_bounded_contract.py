#!/usr/bin/env python3
"""Offline bounded-stop regression/negative controls; never launches MPI."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import yaml

from run import bounded_duration_status, load_module

HERE = Path(__file__).resolve().parent
WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
CASE = WORKSPACE/'results/pending_admission_20261028/pair-r1/full/profile-baseline'
DURATION = load_module(HERE.parent/'full_scale_followup/run.py', 'bounded_contract_reference')


class DurationContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text=(CASE/'aprl.log').read_text()
        cls.record=json.loads((CASE/'aprl.log.launch.json').read_text())

    def check(self,text=None,record=None):
        return bounded_duration_status(self.text if text is None else text,
                                       self.record if record is None else record,20,DURATION)

    def rejected(self,text=None,record=None):
        with self.assertRaises((RuntimeError,ValueError)):
            self.check(text,record)

    def test_observed_18979_post_close_abort(self):
        report=self.check()
        self.assertTrue(report['readiness_complete_marker_observed'])
        self.assertFalse(report['physical_simulation_complete'])
        self.assertEqual(len(report['scheduler_step_cancellations_after_duration_abort']),1)
        self.assertEqual(report['scheduler_step_cancellations_after_duration_abort'][0]['step_id'],'18979.1')

    def test_no_scheduler_cancel_banner_still_valid(self):
        report=self.check('\n'.join(line for line in self.text.splitlines() if 'CANCELLED' not in line))
        self.assertEqual(report['scheduler_step_cancellations_after_duration_abort'],[])

    def test_missing_each_root_close_or_collective_ready_rejected(self):
        for phase in ('trajectory-close-complete','detector-close-complete','ledger-close-complete',
                      'particle-boundary-report-complete','final-resource-report-complete',
                      'duration-exhaustion-readiness-sync-enter','duration-exhaustion-readiness-sync-complete'):
            with self.subTest(phase=phase):
                self.rejected('\n'.join(line for line in self.text.splitlines() if 'phase='+phase not in line))

    def test_timeout_signal_crash_exit_codes_rejected(self):
        for code in (0,124,137,139,143,-9,-15):
            with self.subTest(code=code):
                record=copy.deepcopy(self.record);record['exit_code']=code
                self.rejected(record=record)

    def test_launcher_timeout_even_if_srun_exit_is_one_rejected(self):
        record=copy.deepcopy(self.record);record['seconds']=record['limit_s']
        self.rejected(record=record)

    def test_early_cancel_is_not_misclassified(self):
        cancel=next(line for line in self.text.splitlines() if 'CANCELLED' in line)
        self.rejected(cancel+'\n'+self.text)

    def test_whole_job_or_explained_scheduler_cancel_rejected(self):
        for changed in ('*** JOB 18979','*** STEP 18979.1'):
            if changed.startswith('*** JOB'):
                text=self.text.replace('*** STEP 18979.1',changed)
            else:
                text=self.text.replace('CANCELLED AT','CANCELLED DUE TO TIME LIMIT AT')
            self.rejected(text)

    def test_other_rank_failure_not_hidden_by_expected_messages(self):
        self.rejected(self.text+'\nSimulation failure on MPI rank 3: invalid particle owner\n')

    def test_actual_signal_crash_oom_and_bad_peer_exit_rejected(self):
        for error in ('Segmentation fault','Killed by signal 9','Timed out',
                      'slurmstepd: error: Detected 1 oom_kill event',
                      'srun: error: node5: task 3: Exited with exit code 2'):
            with self.subTest(error=error): self.rejected(self.text+'\n'+error+'\n')

    def test_missing_or_other_mpi_abort_reason_rejected(self):
        self.rejected(self.text.replace('MPI_ABORT was invoked','MPI_ABORT unknown'))
        self.rejected(self.text.replace('with errorcode 1.','with errorcode 2.',1))

    def test_wrong_final_step_rejected(self):
        self.rejected(self.text.replace('phase=final step=20','phase=final step=19'))

    def test_physical_completion_not_promoted(self):
        self.rejected(self.text+'\nConfigured stop reached: all particles passed\n')


class OutputTransactions(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        (self.root/'detectors').mkdir();(self.root/'ledger').mkdir()
        (self.root/'run.yaml').write_text(yaml.safe_dump({'detectors':{
            'field_planes':[{'name':'exit-field'}],'particle_planes':[]},
            'energy_ledger':{'filename':'energy-ledger.h5'}}))
        with h5py.File(self.root/'detectors/exit-field.h5','w') as f:
            g=f.create_group('field_plane');g['complete']=np.uint8(0);g['committed_samples']=np.uint64(2)
            g['time_s']=np.array([0.,1.]);g['fields']=np.zeros((2,1,1))
        with h5py.File(self.root/'ledger/energy-ledger.h5','w') as f:
            g=f.create_group('energy_ledger');g['complete']=np.uint8(0);g['committed_records']=np.uint64(2)
            g['records']=np.array([0.,1.])

    def tearDown(self): self.temp.cleanup()

    def test_closed_complete_false_transaction_accepted(self):
        result=DURATION.inspect_outputs(self.root,False)
        self.assertEqual([item['complete'] for item in result],[0,0])

    def test_uncommitted_or_missing_extent_rejected(self):
        with h5py.File(self.root/'detectors/exit-field.h5','r+') as f:
            f['field_plane/committed_samples'][()]=3
        with self.assertRaises(ValueError): DURATION.inspect_outputs(self.root,False)

    def test_unexpected_complete_true_rejected(self):
        with h5py.File(self.root/'detectors/exit-field.h5','r+') as f:
            f['field_plane/complete'][()]=1
        with self.assertRaises(ValueError): DURATION.inspect_outputs(self.root,False)

    def test_truncated_hdf5_rejected(self):
        path=self.root/'detectors/exit-field.h5'
        with path.open('r+b') as f: f.truncate(64)
        with self.assertRaises(OSError): DURATION.inspect_outputs(self.root,False)

    def test_nonfinite_payload_rejected(self):
        with h5py.File(self.root/'detectors/exit-field.h5','r+') as f:
            f['field_plane/fields'][0,0,0]=np.nan
        with self.assertRaises(ValueError): DURATION.inspect_outputs(self.root,False)


if __name__=='__main__': unittest.main()
