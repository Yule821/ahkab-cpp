#pragma once

#include <Eigen/Dense>

#include "core/circuit.hpp"
#include "results/results.hpp"

namespace ahkab {

struct AcConfig {
  double fstart = 1.0;
  double fstop = 1e6;
  int npts = 100;
  SweepType sweep_type = SweepType::LOGARITHMIC;
};

AcResult run_ac(const Circuit& circuit, const AcConfig& config);

}  // namespace ahkab
