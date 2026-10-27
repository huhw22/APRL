#include "moving_window_remap.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace aprl {
namespace verification {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

std::size_t product(std::size_t a, std::size_t b) {
  if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
    throw std::overflow_error("Moving-window storage size overflow");
  return a * b;
}

YeeComponent& component(EBFieldGrid& fields, unsigned int index) {
  switch (index) {
    case 0: return fields.ex();
    case 1: return fields.ey();
    case 2: return fields.ez();
    case 3: return fields.bx();
    case 4: return fields.by();
    case 5: return fields.bz();
    case 6: return fields.jx();
    case 7: return fields.jy();
    case 8: return fields.jz();
    default: throw std::out_of_range("Moving-window field component index");
  }
}

std::array<std::size_t, 3> shape(const EBGridGeometry& geometry,
                                unsigned int index) {
  const bool electric = index < 3 || index >= 6;
  const unsigned int axis = index % 3;
  const std::size_t cells[3] = {geometry.nx, geometry.ny, geometry.nz};
  std::array<std::size_t, 3> result;
  for (unsigned int d = 0; d < 3; ++d)
    result[d] = cells[d] + (electric ? d != axis : d == axis);
  return result;
}

void validateParameters(const EBGridGeometry& geometry,
                        const EBCPMLParameters& parameters) {
  const std::size_t cells[3] = {geometry.nx, geometry.ny, geometry.nz};
  for (unsigned int d = 0; d < 3; ++d) {
    // Keeping dimensions in this range also makes signed plan arithmetic safe.
    if (cells[d] < 2 || cells[d] > static_cast<std::size_t>(
          std::numeric_limits<int>::max()))
      throw std::invalid_argument("Moving-window dimensions must fit a positive int");
    const std::size_t p = parameters.cells[d];
    if (p != 0 && (p < 2 || p >= cells[d] || p >= cells[d] - p))
      throw std::invalid_argument("Moving-window CPML slabs need a non-PML interior");
  }
  if (!(geometry.dx > 0.0) || !(geometry.dy > 0.0) ||
      !(geometry.dz > 0.0) || !(geometry.dt > 0.0) ||
      !std::isfinite(geometry.dx) || !std::isfinite(geometry.dy) ||
      !std::isfinite(geometry.dz) || !std::isfinite(geometry.dt))
    throw std::invalid_argument("Moving-window geometry must be finite and positive");
  if (!parameters.enabled() || !(parameters.polynomialOrder > 0.0) ||
      !(parameters.targetReflection > 0.0 && parameters.targetReflection < 1.0) ||
      !(parameters.kappaMax >= 1.0) || !(parameters.alphaFraction >= 0.0) ||
      !std::isfinite(parameters.polynomialOrder) ||
      !std::isfinite(parameters.targetReflection) ||
      !std::isfinite(parameters.kappaMax) || !std::isfinite(parameters.alphaFraction))
    throw std::invalid_argument("Moving-window CPML parameters are invalid");
}

struct MoveBlock {
  std::size_t destination, source, values;
};

struct ZeroBlock {
  std::size_t destination, values;
};

void appendMove(std::vector<MoveBlock>& blocks, std::size_t destination,
                std::size_t source, std::size_t values) {
  if (values == 0) return;
  if (!blocks.empty()) {
    MoveBlock& last = blocks.back();
    if (last.destination + last.values == destination &&
        last.source + last.values == source) {
      last.values += values;
      return;
    }
  }
  blocks.push_back({destination, source, values});
}

void appendZero(std::vector<ZeroBlock>& blocks, std::size_t destination,
                std::size_t values) {
  if (values == 0) return;
  if (!blocks.empty() &&
      blocks.back().destination + blocks.back().values == destination) {
    blocks.back().values += values;
    return;
  }
  blocks.push_back({destination, values});
}

} // namespace

struct MovingWindowRemap::Impl {
  struct FieldPlan {
    MoveBlock move;
    std::size_t incomingBeginK, incomingEndK;
  };
  struct HistoryPlan {
    std::vector<MoveBlock> moves;
    std::vector<ZeroBlock> zeros;
  };
  struct Plan {
    int signedCells;
    std::array<FieldPlan, 9> fields;
    std::array<HistoryPlan, 12> histories;
    MovingWindowPlanSummary summary;
  };

  EBGridGeometry geometry;
  std::array<std::array<std::size_t, 3>, 9> shapes;
  std::array<MovingWindowHistoryLayout, 12> layouts;
  std::vector<Plan> plans;
  double setupSeconds = 0.0;

  Impl(const EBGridGeometry& geometryValue, const EBCPMLParameters& parameters,
       const std::vector<int>& supportedShifts) : geometry(geometryValue) {
    validateParameters(geometry, parameters);
    for (unsigned int c = 0; c < 9; ++c) {
      shapes[c] = shape(geometry, c);
      product(product(product(shapes[c][0], shapes[c][1]), shapes[c][2]),
              sizeof(Double));
    }
    const unsigned int pairs[6][2] = {
      {0, 1}, {0, 2}, {1, 2}, {1, 0}, {2, 0}, {2, 1}
    };
    const std::size_t cells[3] = {geometry.nx, geometry.ny, geometry.nz};
    for (unsigned int h = 0; h < 12; ++h) {
      MovingWindowHistoryLayout& layout = layouts[h];
      const bool electric = h >= 6;
      layout.fieldComponent = pairs[h % 6][0] + (electric ? 0 : 3);
      layout.derivativeAxis = pairs[h % 6][1];
      layout.dimensions = shapes[layout.fieldComponent];
      const unsigned int axis = layout.derivativeAxis;
      const std::size_t thickness = parameters.cells[axis];
      layout.activeAxisIndices.reserve(2 * thickness);
      // Production makeProfile uses strict depth>0: integer electric nodes
      // on the slab/interior interface are inactive; magnetic half nodes are
      // active up to the adjacent cell. No conductivity coefficients change.
      for (std::size_t index = 0; index < thickness; ++index)
        layout.activeAxisIndices.push_back(index);
      if (thickness != 0)
        for (std::size_t index = cells[axis] - thickness + (electric ? 1 : 0);
             index < layout.dimensions[axis]; ++index)
          layout.activeAxisIndices.push_back(index);
      const std::size_t cross = axis == 0 ?
        product(layout.dimensions[1], layout.dimensions[2]) :
        (axis == 1 ? product(layout.dimensions[0], layout.dimensions[2]) :
                     product(layout.dimensions[0], layout.dimensions[1]));
      layout.values = product(layout.activeAxisIndices.size(), cross);
      product(layout.values, sizeof(Double));
    }

    plans.reserve(supportedShifts.size());
    for (int shift : supportedShifts) {
      if (shift == 0) continue; // Always available as a strict no-op.
      const std::int64_t signedShift = shift;
      const std::uint64_t magnitude = signedShift < 0 ? -signedShift : signedShift;
      if (magnitude > geometry.nz)
        throw std::invalid_argument("Prepared moving-window shift exceeds window nz");
      if (std::any_of(plans.begin(), plans.end(),
                     [shift](const Plan& plan) { return plan.signedCells == shift; }))
        throw std::invalid_argument("Duplicate prepared moving-window shift");
      plans.push_back(makePlan(shift, static_cast<std::size_t>(magnitude)));
    }
  }

  Plan makePlan(int shift, std::size_t magnitude) const {
    Plan result;
    result.signedCells = shift;
    result.summary.signedCells = shift;
    for (unsigned int c = 0; c < 9; ++c) {
      const auto& dims = shapes[c];
      const std::size_t plane = dims[0] * dims[1];
      FieldPlan& fp = result.fields[c];
      fp.move = {shift < 0 ? magnitude * plane : 0,
                 shift > 0 ? magnitude * plane : 0,
                 (dims[2] - magnitude) * plane};
      fp.incomingBeginK = shift < 0 ? 0 : dims[2] - magnitude;
      fp.incomingEndK = shift < 0 ? magnitude : dims[2];
      result.summary.fieldMovedBytes += fp.move.values * sizeof(Double);
      result.summary.fieldIncomingBytes += magnitude * plane * sizeof(Double);
    }

    for (unsigned int h = 0; h < 12; ++h) {
      const MovingWindowHistoryLayout& layout = layouts[h];
      const unsigned int axis = layout.derivativeAxis;
      HistoryPlan& hp = result.histories[h];
      if (axis != 2) {
        const std::size_t row = layout.dimensions[axis == 0 ? 1 : 0];
        const std::size_t strip = row * layout.dimensions[2];
        const std::size_t incoming = magnitude * row;
        for (std::size_t ordinal = 0;
             ordinal < layout.activeAxisIndices.size(); ++ordinal) {
          const std::size_t base = ordinal * strip;
          appendMove(hp.moves, base + (shift < 0 ? incoming : 0),
                     base + (shift > 0 ? incoming : 0), strip - incoming);
          appendZero(hp.zeros, base + (shift < 0 ? 0 : strip - incoming), incoming);
        }
      } else {
        const std::size_t plane = layout.dimensions[0] * layout.dimensions[1];
        const auto& active = layout.activeAxisIndices;
        for (std::size_t newOrdinal = 0; newOrdinal < active.size(); ++newOrdinal) {
          const std::int64_t oldK = static_cast<std::int64_t>(active[newOrdinal]) + shift;
          auto old = active.end();
          if (oldK >= 0 && static_cast<std::uint64_t>(oldK) < layout.dimensions[2])
            old = std::lower_bound(active.begin(), active.end(),
                                   static_cast<std::size_t>(oldK));
          if (old != active.end() && *old == static_cast<std::size_t>(oldK)) {
            const std::size_t oldOrdinal = static_cast<std::size_t>(old - active.begin());
            appendMove(hp.moves, newOrdinal * plane, oldOrdinal * plane, plane);
          } else {
            appendZero(hp.zeros, newOrdinal * plane, plane);
          }
        }
        // Source ordinal is monotonically above destination for +shift and
        // below for -shift. Reverse the latter copy order before zeroing any
        // invalid slot, so separate blocks cannot overwrite unread sources.
        if (shift < 0) std::reverse(hp.moves.begin(), hp.moves.end());
      }
      for (const MoveBlock& block : hp.moves)
        result.summary.historyMovedBytes += block.values * sizeof(Double);
      for (const ZeroBlock& block : hp.zeros)
        result.summary.historyZeroedBytes += block.values * sizeof(Double);
      result.summary.historyMoveBlocks += hp.moves.size();
      result.summary.historyZeroBlocks += hp.zeros.size();
    }
    return result;
  }

  const Plan& plan(int signedCells) const {
    for (const Plan& candidate : plans)
      if (candidate.signedCells == signedCells) return candidate;
    throw std::invalid_argument("Moving-window shift was not prepared at construction");
  }

  EBConvolutionalPML& validate(EBFieldGrid& fields) const {
    const EBGridGeometry& g = fields.geometry();
    if (g.nx != geometry.nx || g.ny != geometry.ny || g.nz != geometry.nz ||
        g.dx != geometry.dx || g.dy != geometry.dy || g.dz != geometry.dz ||
        g.dt != geometry.dt || g.solver != geometry.solver)
      throw std::invalid_argument("Moving-window field geometry does not match plan");
    for (unsigned int c = 0; c < 9; ++c) {
      const YeeComponent& a = component(fields, c);
      if (a.nx() != shapes[c][0] || a.ny() != shapes[c][1] || a.nz() != shapes[c][2] ||
          a.size() != shapes[c][0] * shapes[c][1] * shapes[c][2])
        throw std::invalid_argument("Moving-window native component shape mismatch");
    }
    auto* cpml = dynamic_cast<EBConvolutionalPML*>(&fields.boundary());
    if (!cpml || cpml->checkpointHistoryCount() != layouts.size())
      throw std::invalid_argument("Moving-window remap requires the existing native CPML");
    for (std::size_t h = 0; h < layouts.size(); ++h)
      if (cpml->checkpointHistory(h).size() != layouts[h].values)
        throw std::invalid_argument("Moving-window CPML history layout mismatch");
    return *cpml;
  }
};

MovingWindowRemap::MovingWindowRemap(
    const EBGridGeometry& geometry, const EBCPMLParameters& parameters,
    const std::vector<int>& supportedShifts)
  : impl_() {
  const Clock::time_point start = Clock::now();
  impl_.reset(new Impl(geometry, parameters, supportedShifts));
  impl_->setupSeconds = elapsed(start);
}

MovingWindowRemap::~MovingWindowRemap() = default;

MovingWindowShiftReport MovingWindowRemap::shift(
    EBFieldGrid& fields, int signedCells, Double newOriginZ,
    MovingWindowIncoming incoming, MovingWindowCurrentPolicy currentPolicy) {
  MovingWindowShiftReport report;
  report.signedCells = signedCells;
  report.newOriginZ = newOriginZ;
  // In particular, do not clear current, validate/rebind boundaries, query a
  // provider, or touch any value on a no-shift call.
  if (signedCells == 0) return report;
  const Clock::time_point allStart = Clock::now();
  if (!std::isfinite(newOriginZ))
    throw std::invalid_argument("Moving-window new origin must be finite");
  if (currentPolicy != MovingWindowCurrentPolicy::Preserve &&
      currentPolicy != MovingWindowCurrentPolicy::Reset)
    throw std::invalid_argument("Unknown moving-window current policy");
  const Impl::Plan& plan = impl_->plan(signedCells);
  EBConvolutionalPML& cpml = impl_->validate(fields);
  const unsigned int count = currentPolicy == MovingWindowCurrentPolicy::Preserve ? 9 : 6;

  const Clock::time_point fieldStart = Clock::now();
  for (unsigned int c = 0; c < count; ++c) {
    Double* data = component(fields, c).data();
    const MoveBlock& move = plan.fields[c].move;
    if (move.values != 0) {
      std::memmove(data + move.destination, data + move.source, move.values * sizeof(Double));
      report.fieldMovedBytes += move.values * sizeof(Double);
      ++report.fieldMoveBlocks;
    }
  }
  if (currentPolicy == MovingWindowCurrentPolicy::Reset)
    for (unsigned int c = 6; c < 9; ++c) {
      YeeComponent& current = component(fields, c);
      std::fill_n(current.data(), current.size(), Double(0.0));
      report.currentResetBytes += current.size() * sizeof(Double);
    }
  report.fieldWallSeconds = elapsed(fieldStart);

  const Clock::time_point historyStart = Clock::now();
  for (std::size_t h = 0; h < impl_->layouts.size(); ++h) {
    Double* data = cpml.checkpointHistory(h).data();
    const Impl::HistoryPlan& hp = plan.histories[h];
    for (const MoveBlock& move : hp.moves)
      std::memmove(data + move.destination, data + move.source, move.values * sizeof(Double));
    // All source values for this term have been read before invalid slots
    // are cleared. Empty disabled-axis histories have no operations.
    for (const ZeroBlock& zero : hp.zeros)
      std::fill_n(data + zero.destination, zero.values, Double(0.0));
  }
  report.historyMovedBytes = plan.summary.historyMovedBytes;
  report.historyZeroedBytes = plan.summary.historyZeroedBytes;
  report.historyDiscardedBytes = report.historyZeroedBytes;
  report.historyMoveBlocks = plan.summary.historyMoveBlocks;
  report.historyZeroBlocks = plan.summary.historyZeroBlocks;
  report.historyWallSeconds = elapsed(historyStart);

  const Clock::time_point incomingStart = Clock::now();
  for (unsigned int c = 0; c < count; ++c) {
    YeeComponent& field = component(fields, c);
    const auto& dims = impl_->shapes[c];
    const Impl::FieldPlan& fp = plan.fields[c];
    const std::size_t incomingValues = (fp.incomingEndK - fp.incomingBeginK) * dims[0] * dims[1];
    report.incomingValues += incomingValues;
    report.fieldIncomingBytes += incomingValues * sizeof(Double);
    if (!incoming.value) {
      std::fill_n(field.data() + fp.incomingBeginK * dims[0] * dims[1],
                  incomingValues, Double(0.0));
    } else {
      for (std::size_t k = fp.incomingBeginK; k < fp.incomingEndK; ++k)
        for (std::size_t j = 0; j < dims[1]; ++j)
          for (std::size_t i = 0; i < dims[0]; ++i)
            field(i, j, k) = incoming.value(incoming.context, c, i, j, k, newOriginZ);
    }
  }
  report.fieldDiscardedBytes = report.fieldIncomingBytes + report.currentResetBytes;
  report.incomingWallSeconds = elapsed(incomingStart);
  report.totalWallSeconds = elapsed(allStart);
  return report;
}

double MovingWindowRemap::setupWallSeconds() const { return impl_->setupSeconds; }

std::size_t MovingWindowRemap::memoryBytes() const {
  std::size_t bytes = sizeof(*this) + sizeof(Impl) + impl_->plans.capacity() * sizeof(Impl::Plan);
  for (const MovingWindowHistoryLayout& layout : impl_->layouts)
    bytes += layout.activeAxisIndices.capacity() * sizeof(std::size_t);
  for (const Impl::Plan& plan : impl_->plans)
    for (const Impl::HistoryPlan& hp : plan.histories) {
      bytes += hp.moves.capacity() * sizeof(MoveBlock);
      bytes += hp.zeros.capacity() * sizeof(ZeroBlock);
    }
  return bytes;
}

std::size_t MovingWindowRemap::preparedShiftCount() const { return impl_->plans.size(); }

MovingWindowPlanSummary MovingWindowRemap::planSummary(int signedCells) const {
  if (signedCells == 0) return MovingWindowPlanSummary();
  return impl_->plan(signedCells).summary;
}

const MovingWindowHistoryLayout& MovingWindowRemap::historyLayout(std::size_t index) const {
  if (index >= impl_->layouts.size())
    throw std::out_of_range("Moving-window history layout index");
  return impl_->layouts[index];
}

} // namespace verification
} // namespace aprl
