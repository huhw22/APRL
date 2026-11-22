"""Run the installed offline tool on the closed native virtual screen."""
import importlib.util
from pathlib import Path
import sys
import numpy as np

case=Path(sys.argv[1]).resolve()
script=case.parents[1]/'source_coupling_extended_20261122_r1/runtime/bin/finite_screen_spectrum.py'
spec=importlib.util.spec_from_file_location('aprl_screen_spectrum',script)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
frequencies=np.unique(np.r_[np.linspace(0,40e12,801),8.19e12,8.19e12/2,8.19e12/3])
module.write_spectra(case/'radiation.h5',case/'spectrum.h5',frequencies)
print('Offline spectrum complete:',len(frequencies),'frequencies',flush=True)
