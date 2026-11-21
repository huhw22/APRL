#ifndef APRL_TRANSVERSE_SOURCE_SHAPE_H
#define APRL_TRANSVERSE_SOURCE_SHAPE_H

#include <array>

namespace aprl { namespace transverse_source_shape {
// Linear node basis convolved with [1/4,1/2,1/4], indexed from floor(x)-1.
inline std::array<double,4> nodeWeights(double f) {
  return {{0.25*(1.-f),0.5-0.25*f,0.25+0.25*f,0.25*f}};
}
} }
#endif
