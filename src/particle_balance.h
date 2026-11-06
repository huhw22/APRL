#ifndef APRL_PARTICLE_BALANCE_H
#define APRL_PARTICLE_BALANCE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// Integer z partitions only. Costs select ownership; they never enter a
// particle, deposition, field, boundary, or energy update.
namespace aprl { namespace particle_balance {

struct LoadStats {
  std::vector<long double> loads;
  long double total=0.,mean=0.,minimum=0.,maximum=0.,imbalance=1.;
};

namespace detail {
inline std::vector<long double> prefix(const std::vector<double>& work) {
  std::vector<long double> sums(work.size()+1,0.);
  for(std::size_t i=0;i<work.size();++i) {
    if(!std::isfinite(work[i]) || work[i]<0.)
      throw std::invalid_argument("Particle-balance costs must be finite and nonnegative");
    sums[i+1]=sums[i]+static_cast<long double>(work[i]);
    if(!std::isfinite(sums[i+1]))
      throw std::overflow_error("Particle-balance total cost overflow");
  }
  return sums;
}
inline bool validCuts(const std::vector<std::size_t>& cuts,std::size_t cells,
    std::size_t ranks,std::size_t minimum,std::size_t maximum) {
  if(cuts.size()!=ranks+1 || cuts.front()!=0 || cuts.back()!=cells)return false;
  for(std::size_t rank=0;rank<ranks;++rank)
    if(cuts[rank]>cuts[rank+1] || cuts[rank+1]-cuts[rank]<minimum ||
       cuts[rank+1]-cuts[rank]>maximum)return false;
  return true;
}
inline std::vector<std::size_t> uniform(std::size_t cells,std::size_t ranks) {
  std::vector<std::size_t> cuts(ranks+1,0);
  const std::size_t base=cells/ranks,remainder=cells%ranks;
  for(std::size_t rank=0;rank<ranks;++rank)
    cuts[rank+1]=cuts[rank]+base+(rank<remainder?1:0);
  return cuts;
}
inline long double maximumLoad(const std::vector<long double>& sums,
    const std::vector<std::size_t>& cuts) {
  long double load=0.;
  for(std::size_t rank=0;rank+1<cuts.size();++rank)
    load=std::max(load,sums[cuts[rank+1]]-sums[cuts[rank]]);
  return load;
}

// With a minimum segment length >1, "take the longest affordable prefix"
// is not a valid feasibility test. For a candidate maximum cost, previous[j]
// means the prefix ending at j has an admissible partition into part-1 slabs.
// Affordable previous endpoints form a sliding suffix: prefix costs are
// nondecreasing. The queue therefore evaluates exact reachability in O(P*N),
// not an approximation based on counts or greedy local cuts.
inline bool feasible(const std::vector<long double>& sums,std::size_t ranks,
    std::size_t minimum,std::size_t maximum,long double limit,
    std::vector<std::size_t>* answer=nullptr) {
  const std::size_t cells=sums.size()-1,missing=std::numeric_limits<std::size_t>::max();
  std::vector<unsigned char> previous(cells+1,0),current(cells+1,0);
  std::vector<std::size_t> queue(cells+1,0),parents;
  if(answer) {
    if(ranks+1>parents.max_size()/(cells+1))
      throw std::overflow_error("Particle-balance reconstruction is too large");
    parents.assign((ranks+1)*(cells+1),missing);
  }
  previous[0]=1;
  for(std::size_t part=1;part<=ranks;++part) {
    std::fill(current.begin(),current.end(),0);
    std::size_t head=0,tail=0;
    const std::size_t first=part*minimum;
    const std::size_t last=std::min(cells-(ranks-part)*minimum,
        maximum>cells/part?cells:part*maximum);
    bool any=false;
    for(std::size_t end=first;end<=last;++end) {
      const std::size_t candidate=end-minimum;
      if(previous[candidate])queue[tail++]=candidate;
      const std::size_t earliest=end>maximum?end-maximum:0;
      while(head<tail && (queue[head]<earliest ||
            sums[end]-sums[queue[head]]>limit))++head;
      if(head<tail) {
        current[end]=1;any=true;
        if(answer)parents[part*(cells+1)+end]=queue[head];
      }
    }
    if(!any)return false;
    previous.swap(current);
  }
  if(!previous[cells])return false;
  if(answer) {
    answer->assign(ranks+1,0);(*answer)[ranks]=cells;
    for(std::size_t part=ranks;part>0;--part)
      (*answer)[part-1]=parents[part*(cells+1)+(*answer)[part]];
  }
  return true;
}
} // namespace detail

inline LoadStats measure(const std::vector<double>& work,
    const std::vector<std::size_t>& cuts) {
  const auto sums=detail::prefix(work);
  if(cuts.size()<2 || !detail::validCuts(cuts,work.size(),cuts.size()-1,1,work.size()))
    throw std::invalid_argument("Particle-balance statistics require a complete nonempty partition");
  LoadStats result;result.total=sums.back();result.mean=result.total/(cuts.size()-1);
  result.minimum=std::numeric_limits<long double>::max();
  for(std::size_t rank=0;rank+1<cuts.size();++rank) {
    const long double load=sums[cuts[rank+1]]-sums[cuts[rank]];
    result.loads.push_back(load);result.minimum=std::min(result.minimum,load);
    result.maximum=std::max(result.maximum,load);
  }
  result.imbalance=result.mean>0.?result.maximum/result.mean:1.;
  return result;
}

// maxCells==0 means no maximum other than the complete global cell count.
// Hysteresis is a fraction of current maximum load: retain a valid current
// partition unless the proposed reduction exceeds this fraction. Ties retain
// the current cuts even when hysteresis==0. A caller supplies its own positive
// field baseline in work; all-zero work is accepted and partitions uniformly.
inline std::vector<std::size_t> plan(const std::vector<double>& work,
    std::size_t ranks,std::size_t minCells=2,std::size_t maxCells=0,
    const std::vector<std::size_t>& currentCuts=std::vector<std::size_t>(),
    double hysteresis=0.) {
  const auto sums=detail::prefix(work);
  const std::size_t cells=work.size();
  if(!ranks || !minCells || cells/minCells<ranks)
    throw std::invalid_argument("Particle-balance ranks/minimum cells are infeasible");
  if(!std::isfinite(hysteresis) || hysteresis<0. || hysteresis>1.)
    throw std::invalid_argument("Particle-balance hysteresis must be in [0,1]");
  const std::size_t maximum=maxCells?std::min(maxCells,cells):cells;
  if(maximum<minCells || maximum<cells/ranks+(cells%ranks?1:0))
    throw std::invalid_argument("Particle-balance maximum cells are infeasible");
  if(!currentCuts.empty() && !detail::validCuts(currentCuts,cells,ranks,minCells,maximum))
    throw std::invalid_argument("Particle-balance current cuts violate partition constraints");
  std::vector<std::size_t> cuts=detail::uniform(cells,ranks);
  if(sums.back()>0.) {
    long double lower=sums.back()/ranks;
    for(double value:work)lower=std::max(lower,static_cast<long double>(value));
    long double upper=detail::maximumLoad(sums,cuts);
    // Upper always corresponds to a feasible explicit uniform partition.
    // Stop once long-double midpoint arithmetic cannot refine the bracket.
    for(unsigned iteration=0;iteration<80;++iteration) {
      const long double middle=lower+(upper-lower)/2;
      if(middle==lower || middle==upper)break;
      if(detail::feasible(sums,ranks,minCells,maximum,middle))upper=middle;
      else lower=middle;
    }
    if(!detail::feasible(sums,ranks,minCells,maximum,upper,&cuts))
      throw std::logic_error("Particle-balance feasible bound lost during reconstruction");
  }
  if(!detail::validCuts(cuts,cells,ranks,minCells,maximum))
    throw std::logic_error("Particle-balance planner produced invalid integer cuts");
  if(!currentCuts.empty()) {
    const long double old=detail::maximumLoad(sums,currentCuts);
    const long double proposed=detail::maximumLoad(sums,cuts);
    const long double slack=32*std::numeric_limits<long double>::epsilon()*std::max(old,proposed);
    if(old-proposed<=static_cast<long double>(hysteresis)*old+slack)return currentCuts;
  }
  return cuts;
}

}} // namespace aprl::particle_balance
#endif
