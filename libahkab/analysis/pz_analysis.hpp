#pragma once

#include "core/circuit.hpp"
#include "results/results.hpp"

namespace ahkab {

struct PzConfig {
  std::string input_source;
  std::string output_node;
};

PzResult run_pz(const Circuit& circuit, const PzConfig& config);

}  // namespace ahkab
