#ifndef APRL_TEST_NATIVE_TWO_MAGNETIC_HALO_H
#define APRL_TEST_NATIVE_TWO_MAGNETIC_HALO_H
// Included inside source_coupling_math in the private test header only.
// Same message count and tags as the original: double Bx/By payloads only.
class NativeHalo {
 bool lower_=false,upper_=false,magnetic_=false;
 std::vector<Double> down_[3],up_[3];
public:
 void exchange(const EBFieldGrid& f,MPI_Comm comm,bool magnetic=false){
   magnetic_=magnetic;int rank=0,size=1;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);lower_=rank>0;upper_=rank+1<size;if(size==1)return;
   const YeeComponent* values[3]={magnetic?&f.bx():&f.ex(),magnetic?&f.by():&f.ey(),magnetic?&f.bz():&f.ez()};
   for(unsigned a=0;a<3;++a){const auto& q=*values[a];const std::size_t plane=q.nx()*q.ny(),width=magnetic&&a<2?2:1,count=width*plane;const bool node=magnetic?a==2:a<2;
     if(count>std::size_t(std::numeric_limits<int>::max())||q.nz()<width+(node?1:0))throw std::overflow_error("Two-plane magnetic collar extent");down_[a].resize(count);up_[a].resize(count);
     MPI_Sendrecv(q.data()+(q.nz()-width-(node?1:0))*plane,int(count),MPI_DOUBLE,upper_?rank+1:MPI_PROC_NULL,11000+2*a,down_[a].data(),int(count),MPI_DOUBLE,lower_?rank-1:MPI_PROC_NULL,11000+2*a,comm,MPI_STATUS_IGNORE);
     MPI_Sendrecv(q.data()+(node?1:0)*plane,int(count),MPI_DOUBLE,lower_?rank-1:MPI_PROC_NULL,11001+2*a,up_[a].data(),int(count),MPI_DOUBLE,upper_?rank+1:MPI_PROC_NULL,11001+2*a,comm,MPI_STATUS_IGNORE);
   }
 }
 bool contains(long long k,std::size_t nz)const{const long long w=magnetic_?2:1;return(k<0&&k>=-w&&lower_)||(k>=(long long)nz&&k<(long long)nz+w&&upper_);}
 bool hasLower()const{return lower_;}bool hasUpper()const{return upper_;}
 Double value(const YeeComponent& v,unsigned a,long long i,long long j,long long k)const{
   if(i<0||j<0||i>=(long long)v.nx()||j>=(long long)v.ny())throw std::out_of_range("Coupling transverse native sample outside grid");if(k>=0&&k<(long long)v.nz())return v(i,j,k);
   const long long width=magnetic_&&a<2?2:1;const std::size_t plane=v.nx()*v.ny(),index=j*v.nx()+i;
   if(k<0&&k>=-width&&lower_)return down_[a][(width+k)*plane+index];
   if(k>=(long long)v.nz()&&k<(long long)v.nz()+width&&upper_)return up_[a][(k-v.nz())*plane+index];
   throw std::out_of_range("Coupling sample outside native MPI collar");
 }
 Double openMagneticValue(const YeeComponent& v,unsigned a,long long i,long long j,long long k)const{
   i=std::max(0LL,std::min(i,(long long)v.nx()-1));j=std::max(0LL,std::min(j,(long long)v.ny()-1));if(k<0&&!lower_)k=0;if(k>=(long long)v.nz()&&!upper_)k=v.nz()-1;return value(v,a,i,j,k);
 }
};
#endif
