#include "finite_radiation_screen.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aprl
{
  namespace
  {
    typedef std::array<double,3> Vec;
    // Same SI constants as the field solver. The analysis kernel can also
    // be built by itself without linking a Maxwell grid or MPI.
    constexpr double c=299792458., epsilon0=8.8541878128e-12;
    constexpr double pi=3.14159265358979323846;
    constexpr double lwScale=1./(4.*pi*epsilon0*c);

    double dot(const Vec& a,const Vec& b)
    {
      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    }

    Vec cross(const Vec& a,const Vec& b)
    {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],
               a[0]*b[1]-a[1]*b[0]}};
    }

    std::size_t product(std::size_t a,std::size_t b)
    {
      if (a && b>std::numeric_limits<std::size_t>::max()/a)
        throw std::length_error("Finite radiation screen dimensions overflow");
      return a*b;
    }

    struct Source
    {
      const RadiationSegment* segment;
      Vec midpoint,du,meanBeta;
      double gamma0,gamma1;
      long double binBase0, binBase1;
    };

    long double baseTime(long double t,double dz,double origin)
    {
      return (t-origin)+static_cast<long double>(std::abs(dz))/c;
    }

    // R = |dz| + transverse^2/(R+|dz|). Keeping the small transverse delay
    // separate avoids subtracting two nearly equal forward path lengths.
    double transverseDelay(double transverse2,double distance,double dz)
    {
      return transverse2/((distance+std::abs(dz))*c);
    }

    long double binArrival(const Vec& source,double x,double y,double z,long double base)
    {
      double dx=x-source[0], dy=y-source[1], dz=z-source[2];
      double transverse2=dx*dx+dy*dy;
      const double distance=std::sqrt(transverse2+dz*dz);
      if (distance==0.) return base;
      return base+transverseDelay(transverse2,distance,dz);
    }

    double inverseDenominator(double gamma,double projection,const Vec& nCrossU)
    {
      // 1/(gamma-n.u) using gamma^2-|u|^2=1 even when beta rounds to one.
      return projection>=0. ? (gamma+projection)/(1.+dot(nCrossU,nCrossU))
                            : 1./(gamma-projection);
    }

    void addTimeBins(double* output,const RadiationScreenSpec& spec,
                     long double start,long double end,const double* area)
    {
      const long double left=start/spec.timeStep, right=end/spec.timeStep;
      if (right<0.L || left>=spec.nt) return;
      // A sub-bin impulse has a known integrated field even when its arrival
      // width is below floating-point resolution at the requested clock origin.
      if (right==left)
        {
          const auto bin=static_cast<std::size_t>(left);
          for (unsigned d=0;d<6;++d) output[6*bin+d]+=area[d]/spec.timeStep;
          return;
        }
      const auto first=static_cast<std::size_t>(std::max(0.L,left));
      // A whole segment inside one bin has overlap/(right-left) exactly one.
      // Keep clipped segments on the general path, including the left edge.
      if (left>=0.L && right<=static_cast<long double>(first+1))
        {
          const double scale=1./spec.timeStep;
          for (unsigned d=0;d<6;++d) output[6*first+d]+=area[d]*scale;
          return;
        }
      const auto past=static_cast<std::size_t>(std::min(
          static_cast<long double>(spec.nt),std::ceil(right)));
      for (std::size_t i=first;i<past;++i)
        {
          const long double overlap=std::min(right,static_cast<long double>(i+1))-
                                    std::max(left,static_cast<long double>(i));
          const double scale=static_cast<double>(overlap/(right-left))/spec.timeStep;
          for (unsigned d=0;d<6;++d) output[6*i+d]+=area[d]*scale;
        }
    }
  }

  FiniteRadiationScreen::FiniteRadiationScreen(const RadiationScreenSpec& spec):spec_(spec)
  {
    if (!spec.nx || !spec.ny || !spec.nt ||
        !std::isfinite(spec.x0) || !std::isfinite(spec.y0) || !std::isfinite(spec.z) ||
        !std::isfinite(spec.dx) || !std::isfinite(spec.dy) || spec.dx<=0. || spec.dy<=0. ||
        !std::isfinite(spec.phaseOrigin) ||
        !std::isfinite(spec.x0+(spec.nx-1)*spec.dx) ||
        !std::isfinite(spec.y0+(spec.ny-1)*spec.dy))
      throw std::invalid_argument("Invalid finite radiation screen geometry, output, or phase origin");
    if (spec.nt && (!std::isfinite(spec.timeStart) || !std::isfinite(spec.timeStep) ||
                    !(spec.timeStep>0.) ||
                    !std::isfinite(spec.timeStart+spec.nt*spec.timeStep)))
      throw std::invalid_argument("Invalid finite radiation screen time bins");
    const std::size_t pixels=product(spec.nx,spec.ny);
    timeField_.resize(product(product(pixels,spec.nt),6));
  }

  std::vector<double> FiniteRadiationScreen::releaseTimeField()
  {
    released_=true;
    std::vector<double> result;
    result.swap(timeField_);
    return result;
  }

  std::size_t FiniteRadiationScreen::storageBytes() const
  {
    return sizeof(double)*(timeField_.capacity()+
                           spec_.frequenciesHz.capacity());
  }

  void FiniteRadiationScreen::accumulate(const std::vector<RadiationSegment>& segments)
  {
    if (released_) throw std::logic_error("Finite radiation screen already finalized");
    std::vector<Source> sources;
    sources.reserve(segments.size());
    for (const auto& segment:segments)
      {
        if (!std::isfinite(segment.t0) || !std::isfinite(segment.t1) ||
            !std::isfinite(segment.charge) || segment.t1<segment.t0)
          throw std::invalid_argument("Invalid lab time or charge in a radiation segment");
        for (unsigned i=0;i<3;++i)
          if (!std::isfinite(segment.x0[i]) || !std::isfinite(segment.x1[i]) ||
              !std::isfinite(segment.u0[i]) || !std::isfinite(segment.u1[i]))
            throw std::invalid_argument("Nonfinite lab state in a radiation segment");
        if (segment.u0==segment.u1 || segment.charge==0.) continue;
        if (segment.t1==segment.t0)
          throw std::invalid_argument("Accelerated radiation segment has zero lab duration");
        Source s;
        s.segment=&segment;
        for (unsigned i=0;i<3;++i)
          {
            s.midpoint[i]=.5*segment.x0[i]+.5*segment.x1[i];
            s.du[i]=segment.u1[i]-segment.u0[i];
          }
        s.gamma0=std::sqrt(1.+dot(segment.u0,segment.u0));
        s.gamma1=std::sqrt(1.+dot(segment.u1,segment.u1));
        if (!std::isfinite(s.gamma0) || !std::isfinite(s.gamma1))
          throw std::invalid_argument("Proper velocity outside the radiation kernel's numeric range");
        for (unsigned i=0;i<3;++i)
          s.meanBeta[i]=(segment.u0[i]+segment.u1[i])/(s.gamma0+s.gamma1);
        s.binBase0=baseTime(segment.t0,spec_.z-segment.x0[2],spec_.timeStart);
        s.binBase1=baseTime(segment.t1,spec_.z-segment.x1[2],spec_.timeStart);
        sources.push_back(s);
      }
    if (sources.empty()) return;

    const auto pixels=spec_.nx*spec_.ny;
    int invalid=0;
    #pragma omp parallel for schedule(static) reduction(|:invalid)
    for (std::size_t p=0;p<pixels;++p)
      {
        const double x=spec_.x0+(p%spec_.nx)*spec_.dx;
        const double y=spec_.y0+(p/spec_.nx)*spec_.dy;
        for (const Source& s:sources)
          {
            const auto& seg=*s.segment;
            Vec n={{x-s.midpoint[0],y-s.midpoint[1],spec_.z-s.midpoint[2]}};
            const double transverse2=n[0]*n[0]+n[1]*n[1], dz=n[2];
            const double distance=std::sqrt(transverse2+dz*dz);
            if (!(distance>0.) || !std::isfinite(distance))
              { invalid|=1; continue; }
            for (auto& component:n) component/=distance;
            const Vec cross0=cross(n,seg.u0), cross1=cross(n,seg.u1);
            const double inv0=inverseDenominator(s.gamma0,dot(n,seg.u0),cross0);
            const double inv1=inverseDenominator(s.gamma1,dot(n,seg.u1),cross1);
            // The difference of u/(gamma-n.u) is needed only transverse to n.
            // Evaluate the small kick directly instead of subtracting its
            // two large field factors when consecutive momenta nearly agree.
            Vec meanMinusN;
            for (unsigned d=0;d<3;++d)
              meanMinusN[d]=s.meanBeta[d]-n[d];
            const double deltaDen=dot(meanMinusN,s.du);
            Vec crossKick=cross(n,s.du);
            const double scale=(seg.charge*lwScale/distance)*inv1;
            for (unsigned d=0;d<3;++d)
              crossKick[d]=(crossKick[d]-cross0[d]*(deltaDen*inv0))*scale;
            const Vec e=cross(n,crossKick);
            const double area[6]={e[0],e[1],e[2],-crossKick[0],-crossKick[1],-crossKick[2]};

            if (spec_.nt)
              {
                long double first=binArrival(seg.x0,x,y,spec_.z,s.binBase0);
                long double last=binArrival(seg.x1,x,y,spec_.z,s.binBase1);
                if (last<first)
                  { invalid|=2; continue; }
                addTimeBins(&timeField_[6*p*spec_.nt],spec_,first,last,area);
              }
          }
      }
    if (invalid&1) throw std::domain_error("Finite radiation screen intersects a source midpoint");
    if (invalid&2) throw std::domain_error("Noncausal lab segment arrival at finite radiation screen");
  }
}
