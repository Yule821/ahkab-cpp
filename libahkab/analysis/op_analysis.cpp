#include "analysis/op_analysis.hpp"

namespace ahkab {

OpResult run_op(const Circuit& circuit, const NrConfig& config) {
  NrResult nr = solve_dc(circuit, config);
  OpResult result;
  result.meta = nr;

  for (const auto& name : circuit.node_names()) {
    if (name == circuit.ground_node()) {
      result.node_voltages["V(" + name + ")"] = 0.0;
    } else {
      int mna_idx = circuit.mna_node_index(name);
      result.node_voltages["V(" + name + ")"] = nr.x(mna_idx);
    }
  }

  int N_eff = circuit.node_count_effective();
  for (size_t i = 0; i < circuit.vsrcs().size(); ++i) {
    result.branch_currents["I(" + circuit.vsrcs()[i].id + ")"] = nr.x(N_eff + i);
  }

  return result;
}

}  // namespace ahkab
