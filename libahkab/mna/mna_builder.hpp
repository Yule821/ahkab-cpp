#pragma once

#include <Eigen/Sparse>

#include "core/circuit.hpp"

namespace ahkab {

struct MnaSystem {
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  int dimension = 0;
  std::map<std::string, int> vsrc_row;
  std::map<std::string, int> ind_row;
};

MnaSystem build_mna_dc(const Circuit& circuit);

}  // namespace ahkab
