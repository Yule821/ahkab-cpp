#pragma once

#include <Eigen/Sparse>

#include "core/types.hpp"

namespace ahkab {
namespace internal {

inline void add_element(Eigen::SparseMatrix<double>& A, int r, int c, double v) {
  A.coeffRef(r, c) += v;
}

inline void add_element_c(Eigen::SparseMatrix<complex_t>& A, int r, int c, complex_t v) {
  A.coeffRef(r, c) += v;
}

}  // namespace internal
}  // namespace ahkab
