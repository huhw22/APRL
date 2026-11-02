#ifndef APRL_TEST_NATIVE_METRIC_HALO_H
#define APRL_TEST_NATIVE_METRIC_HALO_H
// Included inside source_coupling_math in the isolated candidate only.
class NativeHalo {
 bool lower_=false,upper_=false,magnetic_=false;std::vector<Double> down_[3],up_[3];
 unsigned width(unsigned a)const{return magnetic_?(a<2?3:1):(a==2?3:2);}
public:
 void exchange(const EBFieldGrid& f,MPI_Comm comm,bool magnetic=false){magnetic_=magnetic;int rank=0,size=1;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);lower_=rank>0;upper_=rank+1<size;if(size==1)return;
   const YeeComponent* v[3]={magnetic?&f.bx():&f.ex(),magnetic?&f.by():&f.ey(),magnetic?&f.bz():&f.ez()};for(unsigned a=0;a<3;++a){const auto& q=*v[a];const std::size_t plane=q.nx()*q.ny(),w=width(a),count=w*plane;const bool node=magnetic?a==2:a<2;if(count>std::size_t(std::numeric_limits<int>::max())||q.nz()<w+(node?1:0))throw std::runtime_error("Matched collar extent");down_[a].resize(count);up_[a].resize(count);
     MPI_Sendrecv(q.data()+(q.nz()-w-(node?1:0))*plane,int(count),MPI_DOUBLE,upper_?rank+1:MPI_PROC_NULL,11000+2*a,down_[a].data(),int(count),MPI_DOUBLE,lower_?rank-1:MPI_PROC_NULL,11000+2*a,comm,MPI_STATUS_IGNORE);
     MPI_Sendrecv(q.data()+(node?1:0)*plane,int(count),MPI_DOUBLE,lower_?rank-1:MPI_PROC_NULL,11001+2*a,up_[a].data(),int(count),MPI_DOUBLE,upper_?rank+1:MPI_PROC_NULL,11001+2*a,comm,MPI_STATUS_IGNORE);
   }
 }
 bool contains(long long k,std::size_t nz)const{return(k<0&&k>=-3&&lower_)||(k>=(long long)nz&&k<(long long)nz+3&&upper_);}
 bool hasLower()const{return lower_;}bool hasUpper()const{return upper_;}
 Double value(const YeeComponent& v,unsigned a,long long i,long long j,long long k)const{
   if(i<0||j<0||i>=(long long)v.nx()||j>=(long long)v.ny())throw std::out_of_range("Matched transverse collar");if(k>=0&&k<(long long)v.nz())return v(i,j,k);const long long w=width(a);const std::size_t plane=v.nx()*v.ny(),n=j*v.nx()+i;
   if(k<0&&k>=-w&&lower_)return down_[a][(w+k)*plane+n];if(k>=(long long)v.nz()&&k<(long long)v.nz()+w&&upper_)return up_[a][(k-v.nz())*plane+n];throw std::out_of_range("Matched native MPI collar");
 }
 Double openMagneticValue(const YeeComponent& v,unsigned a,long long i,long long j,long long k)const{i=std::max(0LL,std::min(i,(long long)v.nx()-1));j=std::max(0LL,std::min(j,(long long)v.ny()-1));if(k<0&&!lower_)k=0;if(k>=(long long)v.nz()&&!upper_)k=v.nz()-1;return value(v,a,i,j,k);}
};
#endif
