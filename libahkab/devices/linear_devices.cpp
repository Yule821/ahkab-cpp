#include "devices/linear_devices.hpp"

#include <cmath>

#include "devices/stamp_helpers.hpp"

namespace ahkab {

using internal::add_element;
using internal::add_element_c;

void stamp_resistor(const Component& comp, const std::map<std::string, int>& node_idx,
                    Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double g = 1.0 / comp.params.at("R");
  add_element(A, n1, n1, g);
  add_element(A, n1, n2, -g);
  add_element(A, n2, n1, -g);
  add_element(A, n2, n2, g);
}

void stamp_capacitor_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                        Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double gm = 1e-12;
  add_element(A, n1, n1, gm);
  add_element(A, n1, n2, -gm);
  add_element(A, n2, n1, -gm);
  add_element(A, n2, n2, gm);
}

void stamp_capacitor_transient(const Component& comp, const std::map<std::string, int>& node_idx,
                               Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
                               double dt, double vc_prev, double ic_prev,
                               IntegrationMethod method) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double C = comp.params.at("C");
  double geq, ieq;

  if (method == IntegrationMethod::IMPLICIT_EULER) {
    geq = C / dt;
    ieq = geq * vc_prev;
  } else if (method == IntegrationMethod::TRAPEZOIDAL) {
    geq = 2.0 * C / dt;
    ieq = -ic_prev - geq * vc_prev;
  } else {
    double alpha = 0.0;
    switch (method) {
      case IntegrationMethod::GEAR2: alpha = 3.0 / 2.0; break;
      case IntegrationMethod::GEAR3: alpha = 11.0 / 6.0; break;
      case IntegrationMethod::GEAR4: alpha = 25.0 / 12.0; break;
      case IntegrationMethod::GEAR5: alpha = 137.0 / 60.0; break;
      default: alpha = 3.0 / 2.0; break;
    }
    geq = alpha * C / dt;
    ieq = (alpha * C / dt) * vc_prev;
  }

  add_element(A, n1, n1, geq);
  add_element(A, n1, n2, -geq);
  add_element(A, n2, n1, -geq);
  add_element(A, n2, n2, geq);
  b(n1) -= ieq;
  b(n2) += ieq;
}

void stamp_capacitor_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                        Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b,
                        double omega) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double C = comp.params.at("C");
  complex_t Y = complex_t(0.0, omega * C);
  add_element_c(A, n1, n1, Y);
  add_element_c(A, n1, n2, -Y);
  add_element_c(A, n2, n1, -Y);
  add_element_c(A, n2, n2, Y);
}

void stamp_inductor_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                       int ind_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double gm = 1e-12;
  add_element(A, n1, n1, gm);
  add_element(A, n1, n2, -gm);
  add_element(A, n2, n1, -gm);
  add_element(A, n2, n2, gm);
  add_element(A, n1, ind_row, 1.0);
  add_element(A, n2, ind_row, -1.0);
  add_element(A, ind_row, n1, 1.0);
  add_element(A, ind_row, n2, -1.0);
}

void stamp_inductor_transient(const Component& comp, const std::map<std::string, int>& node_idx,
                              int ind_row, Eigen::SparseMatrix<double>& A,
                              Eigen::VectorXd& b, double dt, double il_prev,
                              double vl_prev, IntegrationMethod method) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double L = comp.params.at("L");
  double req, veq;

  if (method == IntegrationMethod::IMPLICIT_EULER) {
    req = L / dt;
    veq = req * il_prev;
  } else if (method == IntegrationMethod::TRAPEZOIDAL) {
    req = 2.0 * L / dt;
    veq = vl_prev + req * il_prev;
  } else {
    double alpha = 0.0;
    switch (method) {
      case IntegrationMethod::GEAR2: alpha = 3.0 / 2.0; break;
      case IntegrationMethod::GEAR3: alpha = 11.0 / 6.0; break;
      case IntegrationMethod::GEAR4: alpha = 25.0 / 12.0; break;
      case IntegrationMethod::GEAR5: alpha = 137.0 / 60.0; break;
      default: alpha = 3.0 / 2.0; break;
    }
    req = alpha * L / dt;
    veq = req * il_prev;
  }

  add_element(A, n1, ind_row, 1.0);
  add_element(A, n2, ind_row, -1.0);
  add_element(A, ind_row, n1, 1.0);
  add_element(A, ind_row, n2, -1.0);
  add_element(A, ind_row, ind_row, -req);
  b(ind_row) -= veq;
}

void stamp_inductor_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                       int ind_row, Eigen::SparseMatrix<complex_t>& A,
                       Eigen::VectorXcd& b, double omega) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  double L = comp.params.at("L");
  complex_t Z = complex_t(0.0, omega * L);
  add_element_c(A, n1, ind_row, 1.0);
  add_element_c(A, n2, ind_row, -1.0);
  add_element_c(A, ind_row, n1, 1.0);
  add_element_c(A, ind_row, n2, -1.0);
  add_element_c(A, ind_row, ind_row, -Z);
}

void stamp_vsource_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                      int vsrc_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n_pos = node_idx.at(comp.node_pins[0]);
  int n_neg = node_idx.at(comp.node_pins[1]);
  add_element(A, n_pos, vsrc_row, 1.0);
  add_element(A, n_neg, vsrc_row, -1.0);
  add_element(A, vsrc_row, n_pos, 1.0);
  add_element(A, vsrc_row, n_neg, -1.0);
  b(vsrc_row) = comp.dc_value;
}

void stamp_vsource_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                      int vsrc_row, Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b) {
  int n_pos = node_idx.at(comp.node_pins[0]);
  int n_neg = node_idx.at(comp.node_pins[1]);
  add_element_c(A, n_pos, vsrc_row, 1.0);
  add_element_c(A, n_neg, vsrc_row, -1.0);
  add_element_c(A, vsrc_row, n_pos, 1.0);
  add_element_c(A, vsrc_row, n_neg, -1.0);
  double phase_rad = comp.ac_phase * M_PI / 180.0;
  b(vsrc_row) = complex_t(comp.ac_magnitude * std::cos(phase_rad),
                           comp.ac_magnitude * std::sin(phase_rad));
}

void stamp_isource_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                      Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n_pos = node_idx.at(comp.node_pins[0]);
  int n_neg = node_idx.at(comp.node_pins[1]);
  b(n_pos) -= comp.dc_value;
  b(n_neg) += comp.dc_value;
}

void stamp_isource_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                      Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b) {
  int n_pos = node_idx.at(comp.node_pins[0]);
  int n_neg = node_idx.at(comp.node_pins[1]);
  double phase_rad = comp.ac_phase * M_PI / 180.0;
  complex_t val(comp.ac_magnitude * std::cos(phase_rad),
                comp.ac_magnitude * std::sin(phase_rad));
  b(n_pos) -= val;
  b(n_neg) += val;
}

void stamp_vcvs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   int vsrc_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int ctrl_n1 = node_idx.at(comp.node_pins[2]);
  int ctrl_n2 = node_idx.at(comp.node_pins[3]);
  double gain = comp.params.at("gain");
  add_element(A, n1, vsrc_row, 1.0);
  add_element(A, n2, vsrc_row, -1.0);
  add_element(A, vsrc_row, n1, 1.0);
  add_element(A, vsrc_row, n2, -1.0);
  add_element(A, vsrc_row, ctrl_n1, -gain);
  add_element(A, vsrc_row, ctrl_n2, gain);
}

void stamp_vcvs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   int vsrc_row, Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int ctrl_n1 = node_idx.at(comp.node_pins[2]);
  int ctrl_n2 = node_idx.at(comp.node_pins[3]);
  double gain = comp.params.at("gain");
  add_element_c(A, n1, vsrc_row, 1.0);
  add_element_c(A, n2, vsrc_row, -1.0);
  add_element_c(A, vsrc_row, n1, 1.0);
  add_element_c(A, vsrc_row, n2, -1.0);
  add_element_c(A, vsrc_row, ctrl_n1, -gain);
  add_element_c(A, vsrc_row, ctrl_n2, gain);
}

void stamp_vccs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int ctrl_n1 = node_idx.at(comp.node_pins[2]);
  int ctrl_n2 = node_idx.at(comp.node_pins[3]);
  double gm = comp.params.at("gain");
  add_element(A, n1, ctrl_n1, gm);
  add_element(A, n1, ctrl_n2, -gm);
  add_element(A, n2, ctrl_n1, -gm);
  add_element(A, n2, ctrl_n2, gm);
}

void stamp_vccs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int ctrl_n1 = node_idx.at(comp.node_pins[2]);
  int ctrl_n2 = node_idx.at(comp.node_pins[3]);
  complex_t gm(comp.params.at("gain"), 0.0);
  add_element_c(A, n1, ctrl_n1, gm);
  add_element_c(A, n1, ctrl_n2, -gm);
  add_element_c(A, n2, ctrl_n1, -gm);
  add_element_c(A, n2, ctrl_n2, gm);
}

void stamp_ccvs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   int vsrc_row, Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int vbranch_row = vsrc_idx.at(comp.node_pins[2]);
  double gain = comp.params.at("gain");
  add_element(A, n1, vsrc_row, 1.0);
  add_element(A, n2, vsrc_row, -1.0);
  add_element(A, vsrc_row, n1, 1.0);
  add_element(A, vsrc_row, n2, -1.0);
  add_element(A, vsrc_row, vbranch_row, -gain);
}

void stamp_ccvs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   int vsrc_row, Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int vbranch_row = vsrc_idx.at(comp.node_pins[2]);
  double gain = comp.params.at("gain");
  add_element_c(A, n1, vsrc_row, 1.0);
  add_element_c(A, n2, vsrc_row, -1.0);
  add_element_c(A, vsrc_row, n1, 1.0);
  add_element_c(A, vsrc_row, n2, -1.0);
  add_element_c(A, vsrc_row, vbranch_row, complex_t(-gain, 0.0));
}

void stamp_cccs_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int vbranch_row = vsrc_idx.at(comp.node_pins[2]);
  double f = comp.params.at("gain");
  add_element(A, n1, vbranch_row, f);
  add_element(A, n2, vbranch_row, -f);
}

void stamp_cccs_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                   const std::map<std::string, int>& vsrc_idx,
                   Eigen::SparseMatrix<complex_t>& A, Eigen::VectorXcd& b) {
  int n1 = node_idx.at(comp.node_pins[0]);
  int n2 = node_idx.at(comp.node_pins[1]);
  int vbranch_row = vsrc_idx.at(comp.node_pins[2]);
  complex_t f(comp.params.at("gain"), 0.0);
  add_element_c(A, n1, vbranch_row, f);
  add_element_c(A, n2, vbranch_row, -f);
}

}  // namespace ahkab
