"""Build observations around the frozen candidate; no production core change."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
def require(v,s):
 if not v:raise RuntimeError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def once(s,a,b):require(s.count(a)==1,'Overlay anchor '+a[:80]);return s.replace(a,b,1)
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--build',required=True,type=Path);p.add_argument('--campaign',choices=['diagnosis','force-test'],default='diagnosis');a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();b=a.build.resolve();r.mkdir(parents=True,exist_ok=False);b.mkdir(parents=True,exist_ok=False);t=Path(__file__).resolve().parent;f=r/'frozen_source';f.mkdir();base=w/'results/cowan_matched_coupling_20261102/build_r3';core=w/'results/cowan_matched_coupling_20261102/build_r2/source/src'
 for q in t.iterdir():
  if q.is_file():shutil.copy2(q,f/q.name)
 for n in ['metric_gather.h','native_metric_halo.h']:shutil.copy2(base/'frozen_source'/n,f/n)
 shutil.copy2(core/'source_coupling_math.h',f/'source_coupling_math.h');shutil.copy2(base/'frozen_source/periodic.cpp',f/'baseline_periodic.cpp');s=(f/'baseline_periodic.cpp').read_text();s=once(s,'const unsigned padding=3,nx=12,ny=10,nz=24;','unsigned padding=3,nx=12,ny=10,nz=24;int forceMode=0;')
 anchor='int main(int argc,char**argv){MPI_Init(&argc,&argv);configureComputeThreads(MPI_THREAD_SINGLE,1);const double dz=1.2543262415359034e-6;EBGridGeometry g(nx+2*padding,ny+2*padding,nz+2*padding,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);'
 replace='''#include "longitudinal_diagnostic.h"
int main(int argc,char**argv){MPI_Init(&argc,&argv);configureComputeThreads(MPI_THREAD_SINGLE,1);
const unsigned refineZ=argc>4?std::stoul(argv[4]):1;const double scaleXY=argc>5?std::stod(argv[5]):1.,phaseZ=argc>6?std::stod(argv[6]):.31,driftUx=argc>7?std::stod(argv[7]):.02;const unsigned extraZ=argc>8?std::stoul(argv[8]):1,numSteps=160*refineZ;const bool spaceCharge=argc>3&&std::string(argv[3])=="spacecharge";nx=std::lround(12*scaleXY);ny=std::lround(10*scaleXY);nz=24*refineZ;forceMode=argc>9?std::stoi(argv[9]):0;if(forceMode>0&&!spaceCharge)throw std::runtime_error("Alternative force only permitted in longitudinal control");const double dz=1.2543262415359034e-6/refineZ;EBGridGeometry g(nx+2*padding,ny+2*padding,nz+2*padding,5e-6/scaleXY,5e-6/scaleXY,dz,dz/SI::c,EBMaxwellSolver::CowanZ);'''
 s=once(s,anchor,replace);lines=s.splitlines(True);i=next(k for k,line in enumerate(lines) if 'const unsigned quiet=argc>2?' in line)
 lines[i]=''' const unsigned quiet=argc>2?std::stoul(argv[2]):1,Px=spaceCharge?nx:6*quiet,Py=spaceCharge?ny:4*quiet,Pz=4*quiet*extraZ;const Double represented=192e5/(2*Px*Py*Pz);std::vector<RelativisticParticleSI> particles;
 for(unsigned iz=0;iz<Pz;++iz)for(unsigned iy=0;iy<Py;++iy)for(unsigned ix=0;ix<Px;++ix)for(int sign:{-1,1}){
  RelativisticParticleSI p;p.id=particles.size();p.charge=sign*SI::elementaryCharge*represented;p.mass=SI::electronMass*represented;p.position[0]=(padding+(ix+.37)*nx/Px)*g.dx;p.position[1]=(padding+(iy+.43)*ny/Py)*g.dy;const double z=(iz+phaseZ)*nz/Pz;p.position[2]=(padding+z+(spaceCharge?sign*.15*refineZ*std::sin(2*pi*z/nz):0.))*g.dz;p.properVelocity[0]=spaceCharge?0.:sign*driftUx;p.properVelocity[1]=(spaceCharge||driftUx==0.)?0.:.003;p.properVelocity[2]=.23;particles.push_back(p);
 }
 if(spaceCharge)longitudinal_diagnostic::initializeSpaceCharge(f,particles);
 longitudinal_diagnostic::Accumulator diagnosis;
''';s=''.join(lines);s=once(s,'step<160','step<numSteps');s=once(s,'maxIterations=std::max(maxIterations,iterations);long double','maxIterations=std::max(maxIterations,iterations);const auto unwrapped=candidate;long double')
 s=once(s,'sumElectricGap+=eg;sumMagneticGap+=bg;sumReconstruction+=reconstruction;previous=end;','sumElectricGap+=eg;sumMagneticGap+=bg;sumReconstruction+=reconstruction;previous=end;const auto diagnostic=longitudinal_diagnostic::observe(f,old,before,unwrapped,usedE,gaussE);diagnosis.add(diagnostic);')
 s=once(s,r'momentum_reconstruction_J,iterations\n',r'momentum_reconstruction_J,iterations,space_gap_J,time_gap_J,Gauss_gap_J,Gauss_charge_L1_relative,continuity_L1_relative,charge_force_adjoint_error_J,frozen_energy_impossible_particles,diagnostic_particles,max_required_transverse_du,max_required_transverse_E_V_per_m,max_required_longitudinal_du\n')
 s=once(s,"<<iterations<<'\\n';","<<iterations;diagnosis.row(rows,diagnostic);")
 s=once(s,r'\"steps\":160',r'\"steps\":"<<numSteps<<"');s=once(s,'MPI_Finalize();','diagnosis.summary(std::cout);MPI_Finalize();return 0;');(f/'periodic_diagnosis.cpp').write_text(s)
 # Test-only force controls: raw original gather, or charge-adjoint Ez.
 # Only 1-D longitudinal cases can select them. Original Maxwell/deposition stay.
 s=s.replace('matchedOrbitElectric(f,&old,Vec(0.),','longitudinal_diagnostic::testOrbitElectric(f,&old,Vec(0.),')
 s=once(s,'std::vector<Vec> guess;','if(forceMode==2)longitudinal_diagnostic::oldChargeForce=longitudinal_diagnostic::nodes(f,old,before).force;std::vector<Vec> guess;')
 s=once(s,'ChargeConservingCurrentDepositor dep(f,Vec(0.));','if(forceMode==2)longitudinal_diagnostic::newChargeForce=longitudinal_diagnostic::nodes(f,copyElectric(f),before).force;ChargeConservingCurrentDepositor dep(f,Vec(0.));')
 s=once(s,'wrapE(f);Double force=0.;','wrapE(f);if(forceMode==2)longitudinal_diagnostic::newChargeForce=longitudinal_diagnostic::nodes(f,copyElectric(f),candidate).force;Double force=0.;')
 s=once(s,'value+=rho*.5L*(m.apply(m.electric[2],e,2,i,j,k)+m.apply(m.electric[2],e,2,i,j,k-1));','value+=rho*(forceMode==0?.5L*(m.apply(m.electric[2],e,2,i,j,k)+m.apply(m.electric[2],e,2,i,j,k-1)):.5L*(e(2,i,j,k)+e(2,i,j,k-1)));')
 s=s.replace('metric.apply(metric.electric[a],E,a,i,j,z)','(forceMode==0?metric.apply(metric.electric[a],E,a,i,j,z):E(a,i,j,z))').replace('metric.apply(metric.magnetic[a],Bplus,a,i,j,z)','(forceMode==0?metric.apply(metric.magnetic[a],Bplus,a,i,j,z):Bplus(a,i,j,z))')
 s=once(s,'metric.apply(metric.electric[a],mean,a,x,y,z)','(forceMode==0?metric.apply(metric.electric[a],mean,a,x,y,z):mean(a,x,y,z))')
 s=once(s,'worstEnergy<1e-7L&&worstWork<1e-8L&&fieldMomentumError<1e-9L','(forceMode==2||(worstEnergy<1e-7L&&worstWork<1e-8L))&&fieldMomentumError<1e-9L')
 (f/'periodic_diagnosis.cpp').write_text(s)
 runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');archive=w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a';common=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_HAS_MOVING_WINDOW=1','-DAPRL_HAS_OPENMP=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(f),'-I'+str(core),'-isystem',str(runtime/'include')];libs=[str(archive),'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-fopenmp','-lfftw3'];commands={};programs={}
 for name in ['periodic_diagnosis','axial_nondispersion']:
  cmd=common+[str(f/(name+'.cpp')),*libs,'-o',str(b/name)];commands[name]=cmd
  with (r/('build-'+name+'.log')).open('w') as out:subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,check=True)
  programs[name]=str(b/name)
 cases=[]
 def case(name,q=4,wave='oblique',rz=1,xy=1,phase=.31,ux=.02,extra=1,force=0):cases.append(dict(name=name,args=[q,wave,rz,xy,phase,ux,extra,force],steps=160*rz,refine_z=rz,scale_xy=xy,wave=wave,quiet=q,force_mode=force,energy_work_identity_expected=force!=2,timeout_s=480))
 case('oblique_q2',q=2);case('oblique_q3',q=3);case('oblique_q4');case('oblique_q4_phase',phase=.61);case('oblique_q4_z2',rz=2);case('oblique_q4_z2_more_particles',rz=2,extra=2);case('oblique_q4_xy_coarser',xy=.5);case('wave_driven_cold_q4',ux=0);case('longitudinal_only_q4',wave='spacecharge',ux=0);case('longitudinal_only_q4_phase',wave='spacecharge',ux=0,phase=.61);case('longitudinal_only_q4_z2',wave='spacecharge',ux=0,rz=2)
 if a.campaign=='force-test':
  cases=[]
  case('longitudinal_raw_current',wave='spacecharge',ux=0,force=1)
  case('longitudinal_charge_force',wave='spacecharge',ux=0,force=2)
  case('longitudinal_matched_z2_more',wave='spacecharge',ux=0,rz=2,extra=2)
  case('longitudinal_raw_current_z2_more',wave='spacecharge',ux=0,rz=2,extra=2,force=1)
  case('longitudinal_charge_force_z2_more',wave='spacecharge',ux=0,rz=2,extra=2,force=2)
  case('longitudinal_raw_dense_z',wave='spacecharge',ux=0,extra=6,force=1)
  case('longitudinal_charge_dense_z',wave='spacecharge',ux=0,extra=6,force=2)
 save(r/'PLAN.json',dict(workspace=str(w),programs=programs,cases=cases,baseline=str(base),production_program=str(w/'build/source-coupling-exit-r3-20261031/aprl'),production_sha256=sha(w/'build/source-coupling-exit-r3-20261031/aprl'),archive_sha256=sha(archive),Cowan_propagator_changed=False,particle_force_changed=any(x["force_mode"]>0 for x in cases),experimental_transverse_redirection_applied=False,all_cases_c_dt_over_dz=1.0,physical_duration_s=160*1.2543262415359034e-6/299792458.,purpose='Independent longitudinal charge/field force, space/time pairing split, and feasibility of transverse compensation only'))
 shutil.copy2(w/'source/aprl/verification/performance/cases/coupled_lab_sbp_audit/launcher.py',f/'launcher.py');shutil.copy2(t/'run.sbatch',r/'run.sbatch');batch=(r/'run.sbatch').read_text().replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r));(r/'run.sbatch').write_text(batch);save(r/'SOURCE_BUILD.json',dict(commands=commands,program_sha256={k:sha(v) for k,v in programs.items()},archive_sha256=sha(archive),original_force_math_sha256=sha(f/'source_coupling_math.h'),matched_gather_sha256=sha(f/'metric_gather.h'),production_rebuilt=False));save(r/'FROZEN.json',{str(p):sha(p) for p in [*f.iterdir(),r/'run.sbatch',r/'PLAN.json',archive,*map(Path,programs.values())] if p.is_file()});print(r/'run.sbatch')
if __name__=='__main__':main()
