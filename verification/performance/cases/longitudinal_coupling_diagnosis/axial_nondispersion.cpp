#define main periodic_diagnosis_main
#include "periodic_diagnosis.cpp"
#undef main
int main(int argc,char** argv){
 MPI_Init(&argc,&argv);configureComputeThreads(MPI_THREAD_SINGLE,1);nx=6;ny=4;
 long double maxTranslation=0.,maxDiscreteReference=0.,referenceTranslation=0.;
 for(unsigned refine:{1u,2u})for(unsigned polarization:{0u,1u})for(int direction:{-1,1}){
  nz=24*refine;const double dz=1.2543262415359034e-6/refine;
  EBGridGeometry g(nx+2*padding,ny+2*padding,nz+2*padding,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);EBFieldGrid f(g);
  std::vector<double> profile(nz);for(unsigned k=0;k<nz;++k)profile[k]=std::sin(.19*k)+.3*std::cos(.71*k)+.15*(k%2?1.:-1.);
  std::vector<long double> er(nz),br(nz);
  for(unsigned k=0;k<nz;++k){er[k]=profile[k];br[k]=(direction==1?profile[(k+1)%nz]:-profile[k])/SI::c;}
  for(unsigned z=padding;z<padding+nz;++z)for(unsigned y=padding;y<padding+ny;++y)for(unsigned x=padding;x<padding+nx;++x){const unsigned k=z-padding;const double b=br[k];if(polarization==0){f.ex()(x,y,z)=profile[k];f.by()(x,y,z)=b;}else{f.ey()(x,y,z)=profile[k];f.bx()(x,y,z)=-b;}}
  wrapE(f);wrapB(f);
  // Independent one-dimensional Maxwell recurrence using the ORIGINAL constants.
  // epsilon0*mu0*c^2 differs from one by 4.35e-14 in the existing program.
  // The Nyquist limit amplifies this fixed baseline difference quadratically.
  const double magneticScale=g.dt/g.dz,electricScale=(g.dt/(SI::mu0*SI::epsilon0))/g.dz;
  for(unsigned step=1;step<=160*refine;++step){
   auto bn=br,en=er;for(unsigned k=0;k<nz;++k)bn[k]-=(long double)magneticScale*(er[(k+1)%nz]-er[k]);
   for(unsigned k=0;k<nz;++k)en[k]-=(long double)electricScale*(bn[k]-bn[(k+nz-1)%nz]);br=bn;er=en;
   f.advanceMagnetic();wrapB(f);f.clearCurrent();f.advanceElectric();wrapE(f);
   for(unsigned z=padding;z<padding+nz;++z){
    const long long k=z-padding;const unsigned expected=((k-direction*(long long)step)%(long long)nz+nz)%nz;
    const double E=polarization==0?f.ex()(padding,padding,z):f.ey()(padding,padding,z),B=polarization==0?f.by()(padding,padding,z):-f.bx()(padding,padding,z);
    const double expectedB=direction==1?profile[(expected+1)%nz]:-profile[expected];
    maxTranslation=std::max(maxTranslation,std::abs((long double)E-profile[expected]));maxTranslation=std::max(maxTranslation,std::abs((long double)SI::c*B-expectedB));
    maxDiscreteReference=std::max(maxDiscreteReference,std::abs((long double)E-er[k]));maxDiscreteReference=std::max(maxDiscreteReference,std::abs((long double)SI::c*(B-br[k])));
    referenceTranslation=std::max(referenceTranslation,std::abs(er[k]-profile[expected]));referenceTranslation=std::max(referenceTranslation,std::abs((long double)SI::c*br[k]-expectedB));
   }
  }
 }
 std::cout<<std::setprecision(21)<<"{\"complete\":"<<(maxDiscreteReference<1e-10L?"true":"false")<<",\"Cowan_propagator_changed\":false,\"c_dt_over_dz\":1,\"directions\":2,\"polarizations\":2,\"longitudinal_refinements\":2,\"max_abs_E_and_cB_translation_error\":"<<maxTranslation<<",\"max_independent_1D_reference_error\":"<<maxDiscreteReference<<",\"original_constants_reference_translation_error\":"<<referenceTranslation<<",\"original_inverse_epsilon_mu_c_squared_minus_one\":"<<(1.L/((long double)SI::epsilon0*SI::mu0*SI::c*SI::c)-1)<<"}\n";
 if(maxDiscreteReference>1e-10L)throw std::runtime_error("Original Cowan does not reproduce independent longitudinal reference");MPI_Finalize();return 0;
}
