#ifndef APRL_TRANSVERSE_REFERENCE_CURRENT_H
#define APRL_TRANSVERSE_REFERENCE_CURRENT_H
#include <cstddef>
#include <memory>
#include <string>
#include "eb_field.h"
namespace aprl {
namespace verification {class BlockedCowanReference;}
// Prescribed equivalent surface current of the finite reference field.
// It depends only on the prepared reference and time, never on total E/B.
// Thus outgoing perturbations retain the original linear Maxwell/CPML update.
class TransverseReferenceCurrent {
public:
  struct Work {long double signedJ=0.,absoluteJ=0.;};
  static void prepare(const std::string& file,
      verification::BlockedCowanReference& reference,const EBGridGeometry& g,
      const std::size_t guards[3],double solverBeta,bool transverseXBinomial3=false);
  TransverseReferenceCurrent(const std::string& file,const EBGridGeometry& reference,
      const std::size_t guards[3],double solverBeta,const std::string& workFile,bool transverseXBinomial3=false);
  ~TransverseReferenceCurrent();
  void apply(EBFieldGrid& fields,std::ptrdiff_t sourceOffset,std::size_t step,
             double translationBeta,bool collectWork);
  Work finish(const EBFieldGrid& fields);
  void close();
  std::size_t memoryBytes() const;
private:
  struct Impl;std::unique_ptr<Impl> impl_;
};
}
#endif
