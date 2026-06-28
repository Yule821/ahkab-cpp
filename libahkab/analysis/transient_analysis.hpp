#pragma once

#include "core/circuit.hpp"
#include "mna/nr_solver.hpp"
#include "results/results.hpp"

namespace ahkab {

struct TranConfig {
  double tstart = 0.0;
  double tstop = 1e-3;
  double tstep = 1e-6;
  double tmax = 1e-4;
  IntegrationMethod method = IntegrationMethod::TRAPEZOIDAL;
  int gear_order = 2;
  double lte_abstol = 1e-6;
  double lte_reltol = 1e-3;
  NrConfig nr;
};

TranResult run_tran(const Circuit& circuit, const TranConfig& config);

}  // namespace ahkab
