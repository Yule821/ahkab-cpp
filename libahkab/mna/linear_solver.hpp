#pragma once

#include <Eigen/Sparse>

namespace ahkab {

Eigen::VectorXd solve_linear(const Eigen::SparseMatrix<double>& A,
                              const Eigen::VectorXd& b);

Eigen::VectorXcd solve_linear_complex(const Eigen::SparseMatrix<std::complex<double>>& A,
                                       const Eigen::VectorXcd& b);

}  // namespace ahkab
