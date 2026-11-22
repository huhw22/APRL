#include "finite_radiation_screen.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace
{
  using aprl::FiniteRadiationScreen;
  using aprl::RadiationScreenSpec;
  using aprl::RadiationSegment;
  typedef std::array<long double,3> Vector;
  const long double pi=std::acos(-1.L), c=299792458.L;
  const long double charge=-1.602176634e-19L, epsilon0=8.8541878128e-12L;
  const long double gamma=1174.171L, beta=std::sqrt(1.L-1.L/(gamma*gamma));
  const long double radius=3.L, zEntrance=.05L, angle=.0008L;
  const long double sourceStart=.2e-9L;

  void require(bool pass,const char* message)
  {
    if (!pass) throw std::runtime_error(message);
  }

  Vector cross(const Vector& a,const Vector& b)
  {
    return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
  }

  long double dot(const Vector& a,const Vector& b)
  {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  }

  Vector position(long double a)
  {
    return {{2.L*radius*std::sin(a/2.L)*std::sin(a/2.L),0.L,
             zEntrance+radius*std::sin(a)}};
  }

  Vector tangent(long double a)
  {
    return {{std::sin(a),0.L,std::cos(a)}};
  }

  RadiationScreenSpec specFor(double z)
  {
    RadiationScreenSpec spec;
    spec.nx=3; spec.ny=2; spec.z=z;
    spec.x0=-55e-6; spec.dx=55e-6; spec.y0=0.; spec.dy=42e-6;
    spec.frequenciesHz={0.,1e12,2e12,300e12};
    spec.phaseOrigin=static_cast<double>(sourceStart+(z-zEntrance)/c);
    spec.timeStart=spec.phaseOrigin-5.e-18;
    spec.timeStep=.25e-18; spec.nt=256;
    return spec;
  }

  std::vector<RadiationSegment> arc(unsigned segments,double delay=0.)
  {
    std::vector<RadiationSegment> data;
    for (unsigned i=0;i<segments;++i)
      {
        const long double a=angle*i/segments, b=angle*(i+1)/segments;
        RadiationSegment seg;
        for (unsigned d=0;d<3;++d)
          {
            seg.x0[d]=position(a)[d]; seg.x1[d]=position(b)[d];
            seg.u0[d]=gamma*beta*tangent(a)[d]; seg.u1[d]=gamma*beta*tangent(b)[d];
          }
        seg.t0=sourceStart+delay+radius*a/(beta*c);
        seg.t1=sourceStart+delay+radius*b/(beta*c);
        seg.charge=charge;
        data.push_back(seg);
      }
    return data;
  }

  struct AnalyticSample
  {
    std::array<long double,6> areaPerRadian;
    long double arrival;
  };

  AnalyticSample analyticAt(long double a,const Vector& observer)
  {
    const Vector r=position(a), t=tangent(a);
    Vector n={{observer[0]-r[0],observer[1]-r[1],observer[2]-r[2]}};
    const long double R=std::sqrt(dot(n,n));
    for(auto& v:n) v/=R;
    Vector velocity={{beta*t[0],beta*t[1],beta*t[2]}};
    Vector gap={{n[0]-velocity[0], n[1]-velocity[1], n[2]-velocity[2]}};
    Vector derivative={{beta*std::cos(a),0.L,-beta*std::sin(a)}};
    long double kappa=1.L-dot(n,velocity);
    const Vector e=cross(n,cross(gap,derivative)), cb=cross(n,e);
    const long double scale=charge/(4*pi*epsilon0*c*R*kappa*kappa);
    AnalyticSample sample;
    for(unsigned d=0;d<3;++d)
      {
        sample.areaPerRadian[d]=scale*e[d];
        sample.areaPerRadian[3+d]=scale*cb[d];
      }
    sample.arrival=sourceStart+(radius*a/beta+R)/c;
    return sample;
  }

  Vector pixelPosition(const RadiationScreenSpec& spec,std::size_t p)
  {
    return {{spec.x0+(p%spec.nx)*spec.dx,spec.y0+(p/spec.nx)*spec.dy,spec.z}};
  }

  // Independent evaluation of the textbook acceleration field in source
  // angle: dT=(1-n.beta)*dt. Composite Simpson integration uses the actual
  // analytic trajectory and no difference of endpoint field factors.
  std::vector<long double> integrate(const RadiationScreenSpec& spec)
  {
    std::vector<long double> expected(spec.nx*spec.ny*6);
    const unsigned intervals=1024;
    for (std::size_t p=0;p<spec.nx*spec.ny;++p)
      for(unsigned i=0;i<=intervals;++i)
        {
          const long double a=angle*i/intervals;
          const AnalyticSample sample=analyticAt(a,pixelPosition(spec,p));
          const long double weight=(i==0 || i==intervals) ? 1.L : (i%2==0 ? 2.L:4.L);
          const long double scale=weight*angle/(3.L*intervals);
          for(unsigned d=0;d<6;++d)
            expected[p*6+d]+=scale*sample.areaPerRadian[d];
        }
    return expected;
  }

  // Integrate all six saved time components and compare with an independent
  // integral of the textbook LW acceleration field along the analytic arc.
  long double areaError(const FiniteRadiationScreen& screen,
                        const std::vector<long double>& ref)
  {
    long double diff=0., total=0.;
    const auto& cfg=screen.spec();
    for(std::size_t p=0;p<cfg.nx*cfg.ny;++p)
      for(unsigned d=0;d<6;++d)
        {
          long double sum=0.;
          for(std::size_t b=0;b<cfg.nt;++b) sum+=screen.timeField()[(p*cfg.nt+b)*6+d];
          const long double exact=ref[p*6+d];
          diff+=std::pow(sum*cfg.timeStep-exact,2);
          total+=static_cast<long double>(exact)*exact;
        }
    return std::sqrt(diff/total);
  }

  // Find each bin edge on the exact arc. This checks the signed field shape
  // on the arrival clock, independently of both the endpoint quadrature and DC.
  long double angleAt(long double arrival,const Vector& observer)
  {
    if(arrival<=analyticAt(0.,observer).arrival) return 0.;
    if(arrival>=analyticAt(angle,observer).arrival) return angle;
    long double lo=0.,hi=angle;
    for(unsigned k=0;k<55;++k)
      {
        const long double mid=(lo+hi)*.5L;
        if(analyticAt(mid,observer).arrival<arrival) lo=mid;
        else hi=mid;
      }
    return (lo+hi)*.5L;
  }

  std::vector<double> analyticBins(const RadiationScreenSpec& spec)
  {
    std::vector<double> bins(spec.nx*spec.ny*spec.nt*6);
    for(std::size_t p=0;p<spec.nx*spec.ny;++p)
      {
        const Vector obs=pixelPosition(spec,p);
        long double a=angleAt(spec.timeStart,obs);
        for(std::size_t i=0;i<spec.nt;++i)
          {
            const long double t=static_cast<long double>(spec.timeStart)+(i+1.L)*spec.timeStep;
            const long double b=angleAt(t,obs);
            if(b>a)
              {
                const unsigned intervals=8;
                const long double h=(b-a)/intervals;
                for(unsigned k=0;k<=intervals;++k)
                  {
                    const auto smp=analyticAt(a+k*h,obs);
                    const unsigned coeff=(k==0 || k==intervals)?1:(k%2==0?2:4);
                    for(unsigned d=0;d<6;++d)
                      bins[(p*spec.nt+i)*6+d]+=coeff*h/(3.L*spec.timeStep)*smp.areaPerRadian[d];
                  }
              }
            a=b;
          }
      }
    return bins;
  }

  void finiteDistance()
  {
    for (double z:{.3524,.5})
      {
        const auto cfg=specFor(z);
        auto timeCfg=cfg;
        timeCfg.frequenciesHz.clear();
        const auto ref=integrate(cfg);
        const auto refTime=analyticBins(cfg);
        long double previous=1.;
        for(unsigned count:{8u,16u,32u,64u})
          {
            FiniteRadiationScreen screen(cfg);
            const auto path=arc(count);
            screen.accumulate(path);
            FiniteRadiationScreen timeOnly(timeCfg);
            timeOnly.accumulate(path);
            require(timeOnly.timeField()==screen.timeField(),
                    "Offline frequency hints changed the saved waveform");
            const long double residual=areaError(screen,ref);
            long double timeDiff=0., norm=0.;
            for(std::size_t k=0;k<refTime.size();++k)
              {
                timeDiff+=std::pow(static_cast<long double>(screen.timeField()[k])-refTime[k],2);
                norm+=static_cast<long double>(refTime[k])*refTime[k];
              }
            const long double timeError=std::sqrt(timeDiff/norm);
            std::cout << "z=" << z << " segments=" << count << " integrated_field_relative_L2="
                      << residual
                      << " time_waveform_relative_L2=" << timeError << '\n';
            require(residual<.32L*previous,"Finite screen arc convergence");
            if(count==64) require(timeError<.002,"Finite screen signed binned waveform mismatch");
            previous=residual;
          }
        require(previous<3e-6L,"Finite-distance acceleration E/cB disagree with analytic trajectory");
      }
  }

  void coherence()
  {
    auto cfg=specFor(.5);
    FiniteRadiationScreen one(cfg), pair(cfg);
    const auto src=arc(16);
    one.accumulate(src);
    pair.accumulate(src);
    pair.accumulate(src);
    auto opposite=src;
    for(auto& seg:opposite) seg.charge=-seg.charge;
    FiniteRadiationScreen cancelled(cfg);
    cancelled.accumulate(src);
    cancelled.accumulate(opposite);
    const auto& x=one.timeField();
    const auto& two=pair.timeField();
    long double scale=0., pairDiff=0., cancelDiff=0.;
    for(std::size_t i=0;i<x.size();++i)
      {
        scale+=static_cast<long double>(x[i])*x[i];
        pairDiff+=std::pow(two[i]-2.L*x[i],2);
        cancelDiff+=std::pow(static_cast<long double>(cancelled.timeField()[i]),2);
      }
    require(scale>0. && std::sqrt(pairDiff/scale)<1e-13L,
            "Coherent equal-charge time fields failed to add");
    require(std::sqrt(cancelDiff/scale)<1e-13L,
            "Opposite-charge time fields failed to cancel");
    // The phase origin belongs to offline processing and cannot change E(t), B(t).
    cfg.phaseOrigin+=.125e-12;
    FiniteRadiationScreen origin(cfg);
    origin.accumulate(src);
    require(origin.timeField()==x,"Offline phase origin changed the physical waveform");
    const std::size_t originalSize=x.size();
    auto freed=one.releaseTimeField();
    require(freed.size()==originalSize && one.timeField().empty(),"Time field handoff failed");
    bool rejected=false;
    try { one.accumulate(src); } catch(const std::logic_error&) { rejected=true; }
    require(rejected,"Post-finalization accumulation was accepted");
  }

  void backwardScreen()
  {
    auto cfg=specFor(.01);
    cfg.phaseOrigin=sourceStart+(zEntrance-cfg.z)/c;
    cfg.timeStart=cfg.phaseOrigin-1e-13;
    cfg.timeStep=2e-13;
    cfg.nt=96;
    FiniteRadiationScreen screen(cfg);
    screen.accumulate(arc(64));
    const long double residual=areaError(screen,integrate(cfg));
    std::cout << "backward_integrated_field_relative_L2=" << residual << '\n';
    require(residual<3e-5L,"Backward radiation polarization or phase disagrees with LW");
  }

  void clippingAndStraight()
  {
    auto cfg=specFor(.3524);
    auto src=arc(1);
    cfg.nx=cfg.ny=1; cfg.x0=cfg.y0=0.;
    FiniteRadiationScreen pulse(cfg);
    pulse.accumulate(src);
    const double dt=cfg.timeStep;
    unsigned middle=0;
    for(unsigned i=0;i+2<cfg.nt;++i)
      if (pulse.timeField()[i*6] && pulse.timeField()[(i+1)*6])
        { middle=i+1; break; }
    require(middle>0,"No interior arrival bin");
    cfg.timeStart+=middle*dt; cfg.nt=1;
    FiniteRadiationScreen slice(cfg);
    slice.accumulate(src);
    const double full=pulse.timeField()[middle*6], clipped=slice.timeField()[0];
    require(std::abs(full-clipped)<1e-5*std::abs(full),"Trimmed arrival bin integral changed");
    src.front().u1=src.front().u0;
    const double interval=src.front().t1-src.front().t0;
    for(unsigned d=0;d<3;++d)
      src.front().x1[d]=src.front().x0[d]+c*interval*src.front().u0[d]/gamma;
    FiniteRadiationScreen straight(cfg);
    straight.accumulate(src);
    for(double x:straight.timeField()) require(x==0.,"Uniform particle generated pulse");
    std::cout << "clipped_bin_relative=" << std::abs(full-clipped)/std::abs(full)
              << " uniform_straight_zero=1\n";
  }
}

int main()
{
  try
    {
      std::cout << std::setprecision(12);
      finiteDistance(); coherence(); clippingAndStraight(); backwardScreen();
      return 0;
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << '\n';
      return 1;
    }
}
