#include "analysis/pz_analysis.hpp"

#include <Eigen/Eigenvalues>

#include "analysis/op_analysis.hpp"

namespace ahkab {

PzResult run_pz(const Circuit& circuit, const PzConfig& config) {
  int N_eff = circuit.node_count_effective();
  auto gnd_idx = circuit.ground_index();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N_eff, N_eff);

  for (const auto& comp : circuit.components()) {
    if (comp.type == ComponentType::RESISTOR) {
      int n1 = circuit.mna_node_index(comp.node_pins[0]);
      int n2 = circuit.mna_node_index(comp.node_pins[1]);
      double g = 1.0 / comp.params.at("R");
      A(n1, n1) -= g;
      A(n1, n2) += g;
      A(n2, n1) += g;
      A(n2, n2) -= g;
    } else if (comp.type == ComponentType::CAPACITOR) {
      int n1 = circuit.mna_node_index(comp.node_pins[0]);
      int n2 = circuit.mna_node_index(comp.node_pins[1]);
      double C = comp.params.at("C");
      A(n1, n1) -= C;
      A(n1, n2) += C;
      A(n2, n1) += C;
      A(n2, n2) -= C;
    }
  }

  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  PzResult result;
  for (int i = 0; i < es.eigenvalues().size(); ++i) {
    if (std::abs(es.eigenvalues()(i).real()) > 1e-15) {
      result.poles.push_back(es.eigenvalues()(i));
    }
  }

  (void)config;
  return result;
}

}  // namespace ahkab
