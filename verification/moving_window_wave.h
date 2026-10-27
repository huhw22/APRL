#ifndef APRL_VERIFICATION_MOVING_WINDOW_WAVE_H
#define APRL_VERIFICATION_MOVING_WINDOW_WAVE_H

#include <string>
#include "moving_window_remap.h"

namespace aprl {
namespace verification {

// Independent finite vacuum-pulse comparison. The large stationary domain
// never supplies fields or histories to the moving domain after t=0.
// The existing Cowan and CPML operators are used unchanged in both runs.
// staticControl instead uses a stationary 296-cell candidate at global origin
// -192, sharing the reference's lower CPML and the same pulse and detector.
// Its different upper boundary geometry is
// a fixed-CPML control, not a strictly single-variable moving-window contrast.
// zCells accepts only the original 8-cell fixture and the production's
// existing 20-cell longitudinal thickness. The 20-cell candidate gains 12
// cells at its upper edge, retaining the same upper CPML/interior interface;
// x6/y0 remains a two-dimensional fixture, not the full production geometry.
// extendedObservation is moving-only: 384 steps and a 768-cell reference at
// origin -384. The candidate geometry, pulse, trajectory and sampling stay
// identical; the default remains the original 176-step/512-cell comparison.
void runMovingWindowWave(const std::string& directory, bool staticControl = false,
                         unsigned zCells = 8,
                         MovingWindowHistoryPolicy historyPolicy = MovingWindowHistoryPolicy::WorldPosition,
                         bool extendedObservation = false);

} // namespace verification
} // namespace aprl

#endif
