#ifndef APRL_Z_PARTITION_H
#define APRL_Z_PARTITION_H

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include <mpi.h>

namespace aprl {

// Integer global cell ownership. cuts[r]..cuts[r+1]-1 belongs to rank r.
// A shared native node is canonically owned by its right-hand cell; the
// physical final node belongs to the final rank. This ownership is for routing
// contributions only: the usual duplicated Yee interface planes remain.
struct ZPartition {
  std::vector<std::size_t> cuts;

  static ZPartition fromCuts(const std::vector<std::size_t>& values) {
    if(values.size()<2 || values.front()!=0)
      throw std::invalid_argument("Invalid integer z partition");
    for(std::size_t r=1;r<values.size();++r)
      if(values[r]<=values[r-1])
        throw std::invalid_argument("Integer z partition needs nonempty ordered slabs");
    ZPartition result;result.cuts=values;return result;
  }

  static ZPartition gather(std::size_t localNz,MPI_Comm communicator) {
    int count=0;
    MPI_Comm_size(communicator,&count);
    if(count<1)throw std::invalid_argument("Empty integer z partition communicator");
    const unsigned long long local=static_cast<unsigned long long>(localNz);
    if(static_cast<std::size_t>(local)!=localNz)
      throw std::overflow_error("Integer z slab exceeds MPI integer representation");
    std::vector<unsigned long long> sizes(static_cast<std::size_t>(count));
    if(MPI_Allgather(&local,1,MPI_UNSIGNED_LONG_LONG,sizes.data(),1,
        MPI_UNSIGNED_LONG_LONG,communicator)!=MPI_SUCCESS)
      throw std::runtime_error("Cannot gather integer z partition");
    std::vector<std::size_t> offsets(static_cast<std::size_t>(count)+1,0);
    for(int r=0;r<count;++r) {
      const auto n=static_cast<std::size_t>(sizes[static_cast<std::size_t>(r)]);
      if(!n || static_cast<unsigned long long>(n)!=sizes[static_cast<std::size_t>(r)] ||
          n>std::numeric_limits<std::size_t>::max()-offsets[static_cast<std::size_t>(r)])
        throw std::invalid_argument("Invalid or overflowing integer z slab size");
      offsets[static_cast<std::size_t>(r)+1]=offsets[static_cast<std::size_t>(r)]+n;
    }
    return fromCuts(offsets);
  }

  std::size_t ranks() const {return cuts.size()-1;}
  std::size_t totalCells() const {return cuts.back();}
  std::size_t offset(int rank) const {return cuts.at(static_cast<std::size_t>(rank));}
  std::size_t cells(int rank) const {
    const auto r=static_cast<std::size_t>(rank);return cuts.at(r+1)-cuts.at(r);
  }
  std::size_t minimumCells() const {
    std::size_t result=std::numeric_limits<std::size_t>::max();
    for(std::size_t r=0;r<ranks();++r)result=std::min(result,cuts[r+1]-cuts[r]);
    return result;
  }
  int ownerCell(std::size_t globalCell) const {
    if(globalCell>=totalCells())throw std::out_of_range("Global z cell outside partition");
    return static_cast<int>(std::upper_bound(cuts.begin(),cuts.end(),globalCell)-cuts.begin()-1);
  }
  int ownerNode(std::size_t globalNode) const {
    if(globalNode>totalCells())throw std::out_of_range("Global z node outside partition");
    return globalNode==totalCells()?static_cast<int>(ranks()-1):ownerCell(globalNode);
  }
};

}
#endif
