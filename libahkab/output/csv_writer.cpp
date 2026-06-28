#include "output/csv_writer.hpp"

#include <cmath>
#include <iomanip>

namespace ahkab {

void write_csv(const TranResult& result, std::ostream& os) {
  os << "time";
  for (const auto& name : result.variable_names) {
    os << "," << name;
  }
  os << "\n";

  for (size_t i = 0; i < result.time.size(); ++i) {
    os << std::scientific << result.time[i];
    for (double v : result.values[i]) {
      os << "," << v;
    }
    os << "\n";
  }
}

void write_table(const OpResult& result, std::ostream& os) {
  os << "DC Operating Point Results\n";
  os << "==========================\n";

  if (result.meta.converged) {
    os << "Converged after " << result.meta.iterations
       << " iterations, residual = " << std::scientific << result.meta.residual
       << "\n\n";
  } else {
    os << "FAILED to converge after " << result.meta.iterations
       << " iterations\n\n";
  }

  os << std::left << std::setw(20) << "Variable"
     << std::right << std::setw(15) << "Value"
     << "\n";
  os << std::string(35, '-') << "\n";

  for (const auto& [name, val] : result.node_voltages) {
    os << std::left << std::setw(20) << name
       << std::right << std::setw(15) << std::scientific << val << "\n";
  }
  for (const auto& [name, val] : result.branch_currents) {
    os << std::left << std::setw(20) << name
       << std::right << std::setw(15) << std::scientific << val << "\n";
  }
}

}  // namespace ahkab
