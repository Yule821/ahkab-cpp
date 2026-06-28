#include "analysis/dc_analysis.hpp"

#include "analysis/op_analysis.hpp"

namespace ahkab {

DcResult run_dc(const Circuit& circuit, const std::string& source_id,
                double start, double stop, double step,
                const NrConfig& config) {
  DcResult result;

  for (const auto& name : circuit.node_names()) {
    result.variable_names.push_back("V(" + name + ")");
  }

  Circuit sweep_circuit = circuit;

  double val = start;
  bool forward = (stop > start);
  while (forward ? (val <= stop + step * 0.5) : (val >= stop - step * 0.5)) {
    for (auto& comp : sweep_circuit.components_mut()) {
      if (comp.id == source_id &&
          (comp.type == ComponentType::VSOURCE ||
           comp.type == ComponentType::ISOURCE)) {
        comp.dc_value = val;
      }
    }

    OpResult op = run_op(sweep_circuit, config);
    result.sweep_values.push_back(val);

    std::vector<double> row;
    for (const auto& name : circuit.node_names()) {
      auto it = op.node_voltages.find("V(" + name + ")");
      if (it != op.node_voltages.end()) {
        row.push_back(it->second);
      } else {
        row.push_back(0.0);
      }
    }
    result.data.push_back(row);

    val += step;
  }

  return result;
}

}  // namespace ahkab
