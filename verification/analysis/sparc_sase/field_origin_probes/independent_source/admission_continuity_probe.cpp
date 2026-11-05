#include "longitudinal_shape.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdlib>
using namespace aprl;
double weight(double z,long long k,bool edge){const auto b=longitudinal_shape::basis(z,edge);return k<b.first||k>=b.first+b.count?0.:b.v[k-b.first];}
double current(double lo,double hi,long long k){
 const double x[4]={-.8611363115940525752,-.3399810435848562648,.3399810435848562648,.8611363115940525752};
 const double w[4]={.3478548451374538574,.6521451548625461426,.6521451548625461426,.3478548451374538574};
 double sum=0.;for(double a=lo;a<hi;){const double b=std::min(hi,std::floor(a)+1.);for(unsigned n=0;n<4;++n)sum+=(b-a)*.5*w[n]*weight((a+b)*.5+(b-a)*.5*x[n],k,true);a=b;}
 return sum;
}
int main(int argc,char**argv){const double beta=argc>1?std::atof(argv[1]):.51479762393903552,begin=argc>2?std::atof(argv[2]):- .2,end=begin+beta;double fullL1=0.,clipL1=0.,retainedL1=0.,maxR=0.,retainedMissingCurrentL1=0.,linearFullL1=0.,linearRetainedL1=0.;
 const auto linearCharge=[](double z,long long k){return std::max(0.,1.-std::abs(z-k));};
 const auto linearCurrent=[](double a,double b,long long k){return std::max(0.,std::min(b,double(k+1))-std::max(a,double(k)));};
 std::cout<<std::setprecision(18)<<"begin_z_cells,end_z_cells,node,rho_begin_times_volume_over_q,rho_end_times_volume_over_q,full_J_dt_area_over_q,clipped_J_dt_area_over_q,full_continuity_residual_over_q,clipped_continuity_residual_over_q\n";
 for(long long k=-8;k<=9;++k){const double d=weight(end,k,false)-weight(begin,k,false),jf=current(begin,end,k),jc=current(0.,end,k),rf=d+jf-current(begin,end,k-1),rc=d+jc-current(0.,end,k-1);fullL1+=std::abs(rf);clipL1+=std::abs(rc);const double ld=linearCharge(end,k)-linearCharge(begin,k);linearFullL1+=std::abs(ld+linearCurrent(begin,end,k)-linearCurrent(begin,end,k-1));if(k>=1){retainedL1+=std::abs(rc);maxR=std::max(maxR,std::abs(rc));retainedMissingCurrentL1+=std::abs(jf-jc);linearRetainedL1+=std::abs(ld+linearCurrent(0.,end,k)-linearCurrent(0.,end,k-1));}std::cout<<begin<<','<<end<<','<<k<<','<<weight(begin,k,false)<<','<<weight(end,k,false)<<','<<jf<<','<<jc<<','<<rf<<','<<rc<<'\n';}
 std::cerr<<std::setprecision(18)<<"full_L1="<<fullL1<<" clipped_L1="<<clipL1<<" retained_nodes_k_ge_1_L1="<<retainedL1<<" max_retained="<<maxR<<" retained_missing_J_dt_area_over_q_L1="<<retainedMissingCurrentL1<<" linear_full_L1="<<linearFullL1<<" linear_clipped_retained_L1="<<linearRetainedL1<<'\n';
 if(fullL1>2e-14||retainedL1<.01)return 1;
}
