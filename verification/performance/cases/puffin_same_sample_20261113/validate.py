"""Checks on closed software fixtures, magnetic probe, and native import."""
import sys
from pathlib import Path

import h5py
import numpy as np

from prepare import read,save,sha


def require(value,message):
    if not value:raise RuntimeError(message)


def smoke(root):
    a=root/'smoke-stock/inputs/3D';b=root/'smoke-matched/inputs/3D'
    checked=[]
    for kind,dataset in [('aperp','aperp'),('electrons','electrons')]:
        afiles=list(a.glob(f'clara_test_{kind}_*.h5'))
        last=max(afiles,key=lambda p:float(h5attr(p,'runInfo','zbarTotal')))
        with h5py.File(last,'r') as f,h5py.File(b/last.name,'r') as g:
            x=f[dataset][:];y=g[dataset][:]
            require(np.isfinite(x).all() and np.isfinite(y).all(),'Nonfinite smoke output')
            require(abs(float(f['runInfo'].attrs['zbarTotal'])*float(f['runInfo'].attrs['Lg'])-.1375)<1e-12,'Incomplete official fixture')
            error=float(np.max(abs(x-y))/max(float(np.max(abs(x))),1e-300))
            require(np.allclose(x,y,rtol=2e-11,atol=2e-13),'Existing Puffin path changed')
            checked.append(dict(dataset=dataset,relative_max_difference=error,shape=list(x.shape)))
    save(root/'SMOKE_VALIDATED.json',dict(passed=True,scope='Official fixed CLARA fixture: original and added-option binaries, unchanged planepole path.',checks=checked))


def h5attr(path,group,key):
    with h5py.File(path,'r') as f:return f[group].attrs[key]


def probe(root):
    info=read(root/'INPUT_VALIDATED.json');scale=info['scale']
    lines=(root/'field_probe.dat').read_text().splitlines()
    s=np.fromstring(next(x for x in lines if x.startswith('SCALE '))[6:],sep=' ')
    names=['pi','qe_C','eps0','me_kg','gain_length_m','cooperation_length_m','npkbar','eta','kappa','lambda_reference_m']
    error=max(abs(s[j]/scale[key]-1) for j,key in enumerate(names))
    require(error<2e-12,'Python/native Fortran scale mismatch')
    a=np.array([np.fromstring(x[6:],sep=' ') for x in lines if x.startswith('FIELD ')])
    z=a[:,0];y=a[:,1];ku=2*np.pi/.028;length=77*.028
    distance=np.minimum(z,length-z);envelope=np.ones(len(z));derivative=np.zeros(len(z));use=distance<.028
    envelope[use]=.5*(1-np.cos(np.pi*distance[use]/.028))
    derivative[use]=.5*np.pi/.028*np.sin(np.pi*distance[use]/.028)*np.where(z[use]<.028,1,-1)
    expected=np.column_stack((np.zeros(len(z)),np.cosh(ku*y)*(envelope*np.sin(ku*z)-derivative/ku*np.cos(ku*z)),
                              np.sinh(ku*y)*envelope*np.cos(ku*z)))
    field_error=float(np.max(abs(a[:,2:]-expected)))
    require(field_error<2e-11,'External magnetic field does not match APRL formula')
    save(root/'FIELD_VALIDATED.json',dict(passed=True,probe_points=len(a),normalized_max_error=field_error,
        coordinate_scale_relative_error=float(error),actual_Fortran_scale=dict(zip(names,map(float,s))),
        scope='Actual compiled custom getBFields vs APRL external field formula, including entrance/exit derivative and off-axis fields.'))


def native(root):
    info=read(root/'INPUT_VALIDATED.json')
    path=root/'import-check/main_electrons_0.h5'
    with h5py.File(root/'particles.h5','r') as f,h5py.File(path,'r') as g:
        a=f['electrons'][:];b=g['electrons'][:]
    require(a.shape==b.shape,'Puffin dropped or added particles')
    ia=np.argsort(a[:,0]);ib=np.argsort(b[:,0])
    require(np.all(np.diff(a[ia,0])>0),'Initial x is not a unique immutable matching key')
    require(np.array_equal(a[ia],b[ib]),'Native import changed actual phase space/weights')
    save(root/'NATIVE_VALIDATED.json',dict(passed=True,records=len(a),all_seven_columns_bitwise=True,
        scope='Closed actual Puffin initial dump, matched by unique x; no assumption of persistent output IDs.',
        input_native_sha256=sha(root/'particles.h5'),physical_charge_C=info['total_charge_C']))


def io(root):
    directory=root/'io-romio321'
    require('IO_PROBE_PASS' in (directory/'probe.log').read_text(),'ROMIO write did not finish')
    with h5py.File(directory/'probe.h5','r') as f:
        data=f['electrons'][:]
    require(data.shape==(1954743,7),'Wrong I/O probe shape')
    require(np.array_equal(data,np.broadcast_to(np.arange(7),data.shape)),'I/O probe data changed')
    ompio=(root/'io-ompio/probe.log').read_text()
    save(root/'IO_VALIDATED.json',dict(passed=True,component='romio321',N=1954743,columns=7,
        all_values_exact=True,ompio_completed='IO_PROBE_PASS' in ompio,
        romio_tail=(directory/'probe.log').read_text().splitlines()[-3:],
        scope='Standalone same-layout HDF5-MPI writer, no particle dynamics or radiation.'))


if __name__=='__main__':
    globals()[sys.argv[1]](Path(sys.argv[2]))
