#pragma once

#include "core/circuit.hpp"
#include "mna/nr_solver.hpp"
#include "results/results.hpp"

namespace ahkab {

OpResult run_op(const Circuit& circuit, const NrConfig& config = NrConfig());

}  // namespace ahkab
