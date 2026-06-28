#include "mna/linear_solver.hpp"

#include <Eigen/SparseLU>

namespace ahkab {

Eigen::VectorXd solve_linear(const Eigen::SparseMatrix<double>& A,
                              const Eigen::VectorXd& b) {
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    return Eigen::VectorXd::Zero(b.size());
  }
  return solver.solve(b);
}

Eigen::VectorXcd solve_linear_complex(const Eigen::SparseMatrix<std::complex<double>>& A,
                                       const Eigen::VectorXcd& b) {
  Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    return Eigen::VectorXcd::Zero(b.size());
  }
  return solver.solve(b);
}

}  // namespace ahkab
