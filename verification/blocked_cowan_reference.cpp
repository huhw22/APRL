#include "blocked_cowan_reference.h"

// With the production option enabled, direct_eb_core provides the common
// implementation. The historical verification-only option can still compile
// that same source without adding FFTW to the default production build.
#ifndef APRL_HAS_MOVING_WINDOW
#include "../src/blocked_cowan_reference.cpp"
#endif
