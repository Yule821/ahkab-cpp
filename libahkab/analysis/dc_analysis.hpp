#pragma once

#include "core/circuit.hpp"
#include "mna/nr_solver.hpp"
#include "results/results.hpp"

namespace ahkab {

DcResult run_dc(const Circuit& circuit, const std::string& source_id,
                double start, double stop, double step,
                const NrConfig& config = NrConfig());

}  // namespace ahkab
