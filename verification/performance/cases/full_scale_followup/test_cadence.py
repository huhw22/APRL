#!/usr/bin/env python3
"""Focused regression for attosecond cadence stored at nanosecond absolute time."""
import unittest

import numpy as np

from run import field_cadence_report


class FieldCadenceTests(unittest.TestCase):
    dt = 1.2543262415359034e-6 / 299792458.0
    gamma = 929.1305492858543
    count = 601
    start = 1.9089993567537292e-8

    def timestamps(self, cadence_scale=1.0):
        offsets = (np.arange(self.count, dtype=np.longdouble) * np.longdouble(self.dt)
                   / np.longdouble(self.gamma) * np.longdouble(cadence_scale))
        return np.asarray(np.longdouble(self.start) + offsets, dtype=np.float64)

    def check(self, times):
        return field_cadence_report(times, self.dt, self.gamma, self.count, -2.00145835e-13)

    def test_accepts_quantized_attosecond_samples(self):
        times = self.timestamps()
        self.assertFalse(np.allclose(np.diff(times), self.dt / self.gamma,
                                    rtol=1e-8, atol=1e-25))
        result = self.check(times)
        self.assertLessEqual(result['maximum_error_fraction_of_bound'], 1)

    def test_rejects_missing_sample(self):
        with self.assertRaises(ValueError):
            self.check(np.delete(self.timestamps(), 300))

    def test_rejects_repeated_sample(self):
        times = self.timestamps()
        times[300] = times[299]
        with self.assertRaises(ValueError):
            self.check(times)

    def test_rejects_sustained_cadence_drift(self):
        with self.assertRaises(ValueError):
            self.check(self.timestamps(1.000001))

    def test_rejects_full_step_jump(self):
        times = self.timestamps()
        times[300:] += self.dt / self.gamma
        with self.assertRaises(ValueError):
            self.check(times)


if __name__ == '__main__':
    unittest.main()
