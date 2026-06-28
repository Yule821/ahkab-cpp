#pragma once

#include <Eigen/Sparse>

#include "core/types.hpp"

namespace ahkab {

void stamp_diode_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                    const DiodeModel& model, Eigen::SparseMatrix<double>& A,
                    Eigen::VectorXd& b, const Eigen::VectorXd& x,
                    double vt);

void stamp_mosfet_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                     const MosfetModel& model, Eigen::SparseMatrix<double>& A,
                     Eigen::VectorXd& b, const Eigen::VectorXd& x,
                     double vt);

void stamp_mosfet_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                     const MosfetModel& model, Eigen::SparseMatrix<complex_t>& A,
                     Eigen::VectorXcd& b, const Eigen::VectorXd& op,
                     double vt, double omega);

void stamp_bjt_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                  const BjtModel& model, Eigen::SparseMatrix<double>& A,
                  Eigen::VectorXd& b, const Eigen::VectorXd& x,
                  double vt);

void stamp_bjt_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                  const BjtModel& model, Eigen::SparseMatrix<complex_t>& A,
                  Eigen::VectorXcd& b, const Eigen::VectorXd& op,
                  double vt, double omega);

}  // namespace ahkab
