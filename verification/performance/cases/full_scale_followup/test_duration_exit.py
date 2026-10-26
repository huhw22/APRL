#!/usr/bin/env python3
"""Offline gates for expected bounded stops; never runs MPI or the simulator."""
from pathlib import Path
import unittest
from run import check_duration_exit, EXHAUSTED

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')


def synthetic(error=None):
    lines = ['[teardown] rank=0 host=nodeA pid=1 elapsed_s=1.0 phase=' + phase
             for phase in ('detector-close-complete', 'ledger-close-complete',
                           'final-resource-report-complete', 'duration-exhaustion-readiness-sync-enter')]
    lines += ['[resource] phase=final step=600 elapsed_s=1.0',
              error or 'Simulation failure on MPI rank 1: ' + EXHAUSTED]
    return '\n'.join(lines)


class DurationExitTests(unittest.TestCase):
    def test_real_accepted_bounded_logs(self):
        root = WORKSPACE / 'results/duration_acceptance/20261025-234447-science-r2'
        for name in ('bounded-1', 'bounded-2'):
            with self.subTest(name=name):
                result = check_duration_exit((root / name / 'aprl.log').read_text(), 1, 200)
                self.assertTrue(result['root_outputs_closed_before_readiness_barrier'])

    def test_old_broken_full_run_still_rejected(self):
        path = WORKSPACE / 'results/full_scale_cowan/20261025-2210-available-nodes/trial/aprl.log'
        with self.assertRaisesRegex(ValueError, 'Missing root output-close'):
            check_duration_exit(path.read_text(), 1)

    def test_truncated_prefix_after_closed_output(self):
        log = synthetic() + '\nSimulation failure on MPI rank srun: error: nodeA: tasks 2: Killed\n'
        result = check_duration_exit(log, 1)
        self.assertEqual(result['truncated_rank_prefixes'], ['Simulation failure on MPI rank '])

    def test_unknown_error_followed_by_expected_same_line_rejected(self):
        log = synthetic('Simulation failure on MPI rank 1: unknown error ' +
                        'Simulation failure on MPI rank 2: ' + EXHAUSTED)
        with self.assertRaisesRegex(ValueError, 'Unexpected complete rank failure'):
            check_duration_exit(log, 1)

    def test_unknown_error_followed_by_expected_next_line_rejected(self):
        log = synthetic('Simulation failure on MPI rank 1: unknown error\n' +
                        'Simulation failure on MPI rank 2: ' + EXHAUSTED)
        with self.assertRaisesRegex(ValueError, 'Unexpected complete rank failure'):
            check_duration_exit(log, 1)

    def test_complete_prefix_partial_message_not_accepted(self):
        log = synthetic() + '\nSimulation failure on MPI rank 2: mesh.duration was exhau'
        with self.assertRaisesRegex(ValueError, 'Unexpected complete rank failure'):
            check_duration_exit(log, 1)

    def test_missing_close_rejected(self):
        log = synthetic().replace('phase=ledger-close-complete', 'phase=ledger-close-enter')
        with self.assertRaisesRegex(ValueError, 'Missing root output-close'):
            check_duration_exit(log, 1)

    def test_timeout_is_not_success(self):
        with self.assertRaises(ValueError):
            check_duration_exit(synthetic(), 124)

    def test_false_final_step_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Missing final complete-step'):
            check_duration_exit(synthetic().replace('step=600', 'step=599'), 1)


if __name__ == '__main__':
    unittest.main()
