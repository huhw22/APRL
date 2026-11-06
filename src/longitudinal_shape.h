#ifndef APRL_LONGITUDINAL_SHAPE_H
#define APRL_LONGITUDINAL_SHAPE_H
#include "eb_field.h"
#include "z_partition.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <mpi.h>
namespace aprl { namespace longitudinal_shape {
// Quintic source and quartic normal-current basis, with the same symmetric
// seven-point continuum-response compensation as the independent 1-D model.
struct Basis {long long first;unsigned count;double v[12];};
enum class Shape {QuinticCompensated,LinearNotch,LinearBinomial3};
inline double notchCoefficient(int offset) {
 const double c[7]={1./64.,-6./64.,15./64.,44./64.,15./64.,-6./64.,1./64.};
 return c[offset+3];
}
inline double binomial3Coefficient(int offset) {
 const double c[3]={.25,.5,.25};
 return c[offset+1];
}
inline Basis basis(double z,bool edge=false,Shape shape=Shape::QuinticCompensated) {
 if(shape==Shape::LinearBinomial3) {
  const double cell=std::floor(z),f=z-cell;Basis r;
  r.first=static_cast<long long>(cell)-1;r.count=edge?3:4;
  std::fill(r.v,r.v+12,0.);
  for(int k=-1;k<=1;++k) {
   const double w=binomial3Coefficient(k);
   if(edge)r.v[k+1]+=w;
   else {r.v[k+1]+=w*(1.-f);r.v[k+2]+=w*f;}
  }
  return r;
 }
 if(shape==Shape::LinearNotch) {
  const double cell=std::floor(z),f=z-cell;Basis r;
  r.first=static_cast<long long>(cell)-3;r.count=edge?7:8;
  std::fill(r.v,r.v+12,0.);
  for(int k=-3;k<=3;++k) {
   const double w=notchCoefficient(k);
   if(edge)r.v[k+3]+=w;
   else {r.v[k+3]+=w*(1.-f);r.v[k+4]+=w*f;}
  }
  return r;
 }
 const unsigned degree=edge?4:5;const double cell=std::floor(z),f=z-cell;
 double n[6]={1.,0.,0.,0.,0.,0.};
 for(unsigned j=1;j<=degree;++j){double saved=0.;for(unsigned k=0;k<j;++k){const double t=n[k]/j;n[k]=saved+(k+1-f)*t;saved=(f+j-k-1)*t;}n[j]=saved;}
 const double p=5.,a=p/6.,b=3*p/40.+p*(p-1)/72.,d=5*p/112.+p*(p-1)/80.+p*(p-1)*(p-2)/1296.;
 const double c[7]={-d/64.,b/16.+3*d/32.,-a/4.-b/4.-15*d/64.,1+a/2.+3*b/8.+5*d/16.,-a/4.-b/4.-15*d/64.,b/16.+3*d/32.,-d/64.};
 Basis r;r.first=static_cast<long long>(cell)-5;r.count=degree+7;std::fill(r.v,r.v+12,0.);
 for(unsigned j=0;j<=degree;++j)for(unsigned k=0;k<7;++k)r.v[j+k]+=n[j]*c[k];return r;
}
inline double nodeWeight(double z,long long node,Shape shape=Shape::QuinticCompensated){const auto b=basis(z,false,shape);return node<b.first||node>=b.first+b.count?0.:b.v[node-b.first];}
// Exact polynomial quadrature on each crossed tensor cell: degree <= 6.
// Callback receives orbit-average current weights, with native Yee indices.
template<class Visit> unsigned weights(const EBGridGeometry& g,const FieldVector<Double>& origin,
 const FieldVector<Double>& begin,const FieldVector<Double>& end,Visit visit,Shape shape=Shape::QuinticCompensated) {
 const double h[3]={g.dx,g.dy,g.dz};double x[3],v[3],cuts[64]={0.,1.};unsigned count=2;
 for(unsigned a=0;a<3;++a){x[a]=(begin[a]-origin[a])/h[a];v[a]=(end[a]-begin[a])/h[a];
  if(v[a]!=0.)for(long long j=static_cast<long long>(std::floor(std::min(x[a],x[a]+v[a])))+1;j<std::max(x[a],x[a]+v[a]);++j){if(count==64)throw std::runtime_error("Quintic orbit exceeds crossing budget");cuts[count++]=(j-x[a])/v[a];}}
 std::sort(cuts,cuts+count);const double q[4]={-.8611363115940525752,-.3399810435848562648,.3399810435848562648,.8611363115940525752},qw[4]={.3478548451374538574,.6521451548625461426,.6521451548625461426,.3478548451374538574};unsigned segments=0;
 for(unsigned s=1;s<count;++s){const double width=cuts[s]-cuts[s-1];if(width<=0.)continue;++segments;
  const unsigned points=(shape==Shape::LinearNotch||shape==Shape::LinearBinomial3)?2:4;
  for(unsigned t=0;t<points;++t){const double node=points==2?(t?0.57735026918962576451:-0.57735026918962576451):q[t],weight=points==2?1.:qw[t];const double u=.5*(cuts[s]+cuts[s-1]+width*node),factor=.5*width*weight,xx=x[0]+u*v[0],yy=x[1]+u*v[1],zz=x[2]+u*v[2];
   long long i=static_cast<long long>(std::floor(xx)),j=static_cast<long long>(std::floor(yy));if(i==static_cast<long long>(g.nx))--i;if(j==static_cast<long long>(g.ny))--j;
   const double fx=xx-i,fy=yy-j;const auto nodes=basis(zz,false,shape),edge=basis(zz,true,shape);
   for(unsigned k=0;k<nodes.count;++k){const double w=factor*nodes.v[k];if(w==0.)continue;for(unsigned b=0;b<2;++b)visit(0,i,j+b,nodes.first+k,w*(b?fy:1.-fy));for(unsigned a=0;a<2;++a)visit(1,i+a,j,nodes.first+k,w*(a?fx:1.-fx));}
   for(unsigned k=0;k<edge.count;++k){const double w=factor*edge.v[k];if(w==0.)continue;for(unsigned b=0;b<2;++b)for(unsigned a=0;a<2;++a)visit(2,i+a,j+b,edge.first+k,w*(a?fx:1.-fx)*(b?fy:1.-fy));}
  }
 }return segments;
}
// z ghost source contributions are reduced before the existing shared-plane
// current sum. Shared node planes themselves retain the original ownership.
class Current {
 EBFieldGrid& f_;MPI_Comm comm_;int rank_,size_;unsigned depth;Shape shape_;
 std::vector<double> lower_[3],upper_[3],receive_;
 ZPartition partition_;bool thin_;
 void completeThin() {
  const auto displacement=[](const std::vector<int>& counts,std::vector<int>& offsets)->std::size_t {
   std::size_t total=0;offsets.resize(counts.size());
   for(std::size_t p=0;p<counts.size();++p){if(counts[p]<0||total>static_cast<std::size_t>(std::numeric_limits<int>::max())-static_cast<std::size_t>(counts[p]))throw std::overflow_error("Current plane routing exceeds MPI count");offsets[p]=static_cast<int>(total);total+=counts[p];}return total;
  };
  if(partition_.totalCells()>static_cast<std::size_t>(std::numeric_limits<long long>::max()))
   throw std::overflow_error("Current global z indices exceed signed routing range");
  for(unsigned a=0;a<3;++a) {
   auto& v=component(a);const std::size_t plane=v.nx()*v.ny();const bool node=a<2;
   const long long last=static_cast<long long>(partition_.totalCells())-(node?0:1);
   std::vector<std::vector<unsigned long long> > indices(static_cast<std::size_t>(size_));
   std::vector<std::vector<std::size_t> > sources(static_cast<std::size_t>(size_));
   for(unsigned side=0;side<2;++side)for(unsigned k=0;k<depth;++k) {
    const long long global=static_cast<long long>(partition_.offset(rank_))+
      (side?static_cast<long long>(v.nz()+k):static_cast<long long>(k)-depth);
    if(global<0||global>last)continue;
    const int owner=node?partition_.ownerNode(static_cast<std::size_t>(global)):
      partition_.ownerCell(static_cast<std::size_t>(global));
    indices[owner].push_back(static_cast<unsigned long long>(global));
    sources[owner].push_back(static_cast<std::size_t>(side)*depth+k);
   }
   std::vector<int> sendCounts(static_cast<std::size_t>(size_)),receiveCounts(static_cast<std::size_t>(size_));
   for(int p=0;p<size_;++p)sendCounts[p]=static_cast<int>(indices[p].size());
   if(MPI_Alltoall(sendCounts.data(),1,MPI_INT,receiveCounts.data(),1,MPI_INT,comm_)!=MPI_SUCCESS)
    throw std::runtime_error("Cannot exchange current collar counts");
   std::vector<int> sendOffsets,receiveOffsets;
   const auto sent=displacement(sendCounts,sendOffsets),received=displacement(receiveCounts,receiveOffsets);
   std::vector<unsigned long long> sendIndices(sent),receiveIndices(received);
   std::vector<std::size_t> source(sent);
   for(int p=0;p<size_;++p){std::copy(indices[p].begin(),indices[p].end(),sendIndices.begin()+sendOffsets[p]);std::copy(sources[p].begin(),sources[p].end(),source.begin()+sendOffsets[p]);}
   if(MPI_Alltoallv(sendIndices.data(),sendCounts.data(),sendOffsets.data(),MPI_UNSIGNED_LONG_LONG,
       receiveIndices.data(),receiveCounts.data(),receiveOffsets.data(),MPI_UNSIGNED_LONG_LONG,comm_)!=MPI_SUCCESS)
    throw std::runtime_error("Cannot route current collar indices");
   std::vector<int> valueSendCounts(static_cast<std::size_t>(size_)),valueReceiveCounts(static_cast<std::size_t>(size_));
   for(int p=0;p<size_;++p){if(static_cast<std::size_t>(sendCounts[p])>static_cast<std::size_t>(std::numeric_limits<int>::max())/plane||static_cast<std::size_t>(receiveCounts[p])>static_cast<std::size_t>(std::numeric_limits<int>::max())/plane)throw std::overflow_error("Current collar reply exceeds MPI count");valueSendCounts[p]=static_cast<int>(sendCounts[p]*plane);valueReceiveCounts[p]=static_cast<int>(receiveCounts[p]*plane);}
   std::vector<int> valueSendOffsets,valueReceiveOffsets;
   const auto valueSent=displacement(valueSendCounts,valueSendOffsets),valueReceived=displacement(valueReceiveCounts,valueReceiveOffsets);
   std::vector<double> sendValues(valueSent),receiveValues(valueReceived);
   for(std::size_t n=0;n<sent;++n){const auto slot=source[n];const auto& values=slot<depth?lower_[a]:upper_[a];const auto k=slot%depth;std::copy(values.begin()+k*plane,values.begin()+(k+1)*plane,sendValues.begin()+n*plane);}
   if(MPI_Alltoallv(sendValues.data(),valueSendCounts.data(),valueSendOffsets.data(),MPI_DOUBLE,
       receiveValues.data(),valueReceiveCounts.data(),valueReceiveOffsets.data(),MPI_DOUBLE,comm_)!=MPI_SUCCESS)
    throw std::runtime_error("Cannot route current collar planes");
   // Only ghost contributions move. Existing local copies of a shared node
   // stay untouched, and the original tangential shared-plane sum follows.
   for(std::size_t n=0;n<received;++n){const auto global=static_cast<std::size_t>(receiveIndices[n]);if(global<partition_.offset(rank_)||global-partition_.offset(rank_)>=v.nz())throw std::logic_error("Current collar reached a non-owning rank");const auto k=global-partition_.offset(rank_);for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)v(i,j,k)+=receiveValues[n*plane+j*v.nx()+i];}
  }
 }
public:
 Current(EBFieldGrid& f,MPI_Comm c,Shape shape=Shape::QuinticCompensated,unsigned haloDepth=8):f_(f),comm_(c),depth(haloDepth),shape_(shape),partition_(ZPartition::gather(f.geometry().nz,c)),thin_(partition_.minimumCells()<depth){if(!depth||(depth<8&&(shape_!=Shape::LinearBinomial3||depth<2)))throw std::invalid_argument("Compact current collar requires binomial3 and at least two layers");MPI_Comm_rank(c,&rank_);MPI_Comm_size(c,&size_);for(unsigned a=0;a<3;++a){auto& v=component(a);const auto n=depth*v.nx()*v.ny();if(n>std::numeric_limits<int>::max())throw std::overflow_error("Quintic current halo exceeds MPI count");lower_[a].resize(n);upper_[a].resize(n);}}
 YeeComponent& component(unsigned a){return a==0?f_.jx():(a==1?f_.jy():f_.jz());}
 void clear(){for(unsigned a=0;a<3;++a){std::fill(lower_[a].begin(),lower_[a].end(),0.);std::fill(upper_[a].begin(),upper_[a].end(),0.);}}
 unsigned deposit(const FieldVector<Double>& origin,const FieldVector<Double>& begin,const FieldVector<Double>& end,double charge){const auto& g=f_.geometry();
  if(depth<8){const double distance=std::abs(end[2]-begin[2]),tolerance=128*std::numeric_limits<double>::epsilon()*std::max(g.dz,std::max(std::abs(origin[2]),std::max(std::abs(begin[2]),std::abs(end[2]))));if(!std::isfinite(distance)||distance>g.dz+tolerance)throw std::out_of_range("Compact binomial3 current requires at most one-cell longitudinal motion");}
  double scale[3];for(unsigned a=0;a<3;++a)scale[a]=charge*(end[a]-begin[a])/(g.dx*g.dy*g.dz*g.dt);
  return weights(g,origin,begin,end,[&](unsigned a,long long i,long long j,long long k,double w){if(scale[a]==0.)return;auto& v=component(a);if(i<0||j<0||i>=static_cast<long long>(v.nx())||j>=static_cast<long long>(v.ny()))throw std::out_of_range("Quintic transverse current support outside grid");const auto plane=v.nx()*v.ny();const double value=scale[a]*w;if(k>=0&&k<static_cast<long long>(v.nz()))v(i,j,k)+=value;else if(k<0&&rank_>0){if(k<-(long long)depth)throw std::out_of_range("Quintic lower current collar");lower_[a][(depth+k)*plane+j*v.nx()+i]+=value;}else if(k>=static_cast<long long>(v.nz())&&rank_+1<size_){if(k>=static_cast<long long>(v.nz()+depth))throw std::out_of_range("Quintic upper current collar");upper_[a][(k-v.nz())*plane+j*v.nx()+i]+=value;}},shape_);
 }
 void complete(){if(thin_){completeThin();return;}for(unsigned a=0;a<3;++a){auto& v=component(a);const auto plane=v.nx()*v.ny(),n=lower_[a].size();receive_.assign(n,0.);const int lo=rank_?rank_-1:MPI_PROC_NULL,hi=rank_+1<size_?rank_+1:MPI_PROC_NULL;
   MPI_Sendrecv(lower_[a].data(),n,MPI_DOUBLE,lo,22100+a,receive_.data(),n,MPI_DOUBLE,hi,22100+a,comm_,MPI_STATUS_IGNORE);
   if(hi!=MPI_PROC_NULL)for(unsigned k=0;k<depth;++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)v(i,j,f_.geometry().nz-depth+k)+=receive_[k*plane+j*v.nx()+i];
   std::fill(receive_.begin(),receive_.end(),0.);MPI_Sendrecv(upper_[a].data(),n,MPI_DOUBLE,hi,22110+a,receive_.data(),n,MPI_DOUBLE,lo,22110+a,comm_,MPI_STATUS_IGNORE);
   if(lo!=MPI_PROC_NULL)for(unsigned k=0;k<depth;++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)v(i,j,k+(a<2?1:0))+=receive_[k*plane+j*v.nx()+i];
  }}
};
}}
#endif
