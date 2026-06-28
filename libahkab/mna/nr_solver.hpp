#pragma once

#include <Eigen/Dense>

#include "core/circuit.hpp"
#include "core/types.hpp"

namespace ahkab {

struct NrResult {
  Eigen::VectorXd x;
  int iterations = 0;
  double residual = 0.0;
  bool converged = false;
};

struct NrConfig {
  double abstol = 1e-12;
  double reltol = 1e-6;
  int max_iter = 300;
  double gmin = 1e-12;
  int gmin_steps = 10;
};

NrResult solve_dc(const Circuit& circuit, const NrConfig& config);

}  // namespace ahkab
