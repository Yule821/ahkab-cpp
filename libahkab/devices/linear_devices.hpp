#pragma once

#include <Eigen/Sparse>

#include "core/types.hpp"

namespace ahkab {

void stamp_resistor(const Component& comp, const std::map<std::string, int>& node_idx,
                    Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_capacitor_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                        Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_capacitor_transient(const Component& comp, const std::map<std::string, int>& node_idx,
                               Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
                               double dt, double vc_prev, double ic_prev,
                               IntegrationMethod method);

void stamp_capacitor_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                        Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b,
                        double omega);

void stamp_inductor_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                       int ind_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_inductor_transient(const Component& comp, const std::map<std::string, int>& node_idx,
                              int ind_row, Eigen::SparseMatrix<double>& A,
                              Eigen::VectorXd& b, double dt, double il_prev,
                              double vl_prev, IntegrationMethod method);

void stamp_inductor_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                       int ind_row, Eigen::SparseMatrix<complex_t>& A,
                       Eigen::VectorXcd& b, double omega);

void stamp_vsource_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                      int vsrc_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_vsource_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                      int vsrc_row, Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b);

void stamp_isource_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                      Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_isource_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                      Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b);

void stamp_vcvs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   int vsrc_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_vcvs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   int vsrc_row, Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b);

void stamp_vccs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_vccs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b);

void stamp_ccvs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   int vsrc_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_ccvs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   int vsrc_row, Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b);

void stamp_cccs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);

void stamp_cccs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b);

}  // namespace ahkab
