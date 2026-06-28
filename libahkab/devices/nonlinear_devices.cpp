#include "devices/nonlinear_devices.hpp"

#include <algorithm>
#include <cmath>

#include "devices/stamp_helpers.hpp"

namespace ahkab {

using internal::add_element;
using internal::add_element_c;

void stamp_diode_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                    const DiodeModel& model, Eigen::SparseMatrix<double>& A,
                    Eigen::VectorXd& b, const Eigen::VectorXd& x,
                    double vt) {
  int na = node_idx.at(comp.node_pins[0]);
  int nc = node_idx.at(comp.node_pins[1]);
  double vd = x(na) - x(nc);

  double nvt = model.n * vt;
  double vcrit = nvt * std::log(nvt / (std::sqrt(2.0) * model.is));
  double id, gd;

  if (vd >= vcrit + 10.0 * nvt) {
    double evd = std::exp(vcrit / nvt);
    double evd_crit = std::exp((vcrit - vcrit) / nvt);
    id = model.is * evd * (1.0 + (vd - vcrit) / nvt);
    gd = model.is * evd / nvt;
  } else if (vd < -5.0 * nvt) {
    id = -model.is;
    gd = 0.0;
  } else {
    double evd = std::exp(vd / nvt);
    id = model.is * (evd - 1.0);
    gd = model.is * evd / nvt;
  }

  double ieq = id - gd * vd;
  add_element(A, na, na, gd);
  add_element(A, na, nc, -gd);
  add_element(A, nc, na, -gd);
  add_element(A, nc, nc, gd);
  b(na) -= ieq;
  b(nc) += ieq;
}

void stamp_mosfet_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                     const MosfetModel& model, Eigen::SparseMatrix<double>& A,
                     Eigen::VectorXd& b, const Eigen::VectorXd& x,
                     double vt) {
  int nd = node_idx.at(comp.node_pins[0]);
  int ng = node_idx.at(comp.node_pins[1]);
  int ns = node_idx.at(comp.node_pins[2]);
  int nb = node_idx.at(comp.node_pins[3]);

  double vgs = x(ng) - x(ns);
  double vds = x(nd) - x(ns);
  double vbs = x(nb) - x(ns);

  double vth = model.vto + model.gamma * (std::sqrt(std::max(0.0, model.phi - vbs)) -
                                          std::sqrt(model.phi));
  double vgst = vgs - vth;

  double ids, gm, gds, gmbs;
  if (vgst <= 0.0) {
    ids = 0.0;
    gm = 0.0;
    gds = 0.0;
    gmbs = 0.0;
  } else if (vds <= vgst) {
    double k = model.kp;
    ids = k * (vgst * vds - 0.5 * vds * vds) * (1.0 + model.lambda * vds);
    gm = k * vds * (1.0 + model.lambda * vds);
    gds = k * (vgst - vds) * (1.0 + model.lambda * vds) + k * (vgst * vds - 0.5 * vds * vds) * model.lambda;
    gmbs = -model.gamma * 0.5 / std::sqrt(std::max(1e-12, model.phi - vbs)) * gm;
  } else {
    double k = model.kp;
    double vdsat = vgst;
    ids = 0.5 * k * vgst * vgst * (1.0 + model.lambda * vds);
    gm = k * vgst * (1.0 + model.lambda * vds);
    gds = 0.5 * k * vgst * vgst * model.lambda;
    gmbs = -model.gamma * 0.5 / std::sqrt(std::max(1e-12, model.phi - vbs)) * gm;
  }

  double ids_eq = ids - gm * vgs - gds * vds - gmbs * vbs;

  add_element(A, nd, nd, gds);
  add_element(A, nd, ns, -gds - gm - gmbs);
  add_element(A, nd, ng, gm);
  add_element(A, nd, nb, gmbs);
  add_element(A, ns, nd, -gds);
  add_element(A, ns, ns, gds + gm + gmbs);
  add_element(A, ns, ng, -gm);
  add_element(A, ns, nb, -gmbs);

  b(nd) -= ids_eq;
  b(ns) += ids_eq;
}

void stamp_mosfet_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                     const MosfetModel& model, Eigen::SparseMatrix<complex_t>& A,
                     Eigen::VectorXcd& b, const Eigen::VectorXd& op,
                     double vt, double omega) {
  int nd = node_idx.at(comp.node_pins[0]);
  int ng = node_idx.at(comp.node_pins[1]);
  int ns = node_idx.at(comp.node_pins[2]);
  int nb = node_idx.at(comp.node_pins[3]);

  double vgs = op(ng) - op(ns);
  double vds = op(nd) - op(ns);
  double vbs = op(nb) - op(ns);

  double vth = model.vto + model.gamma * (std::sqrt(std::max(0.0, model.phi - vbs)) -
                                          std::sqrt(model.phi));
  double vgst = vgs - vth;

  double gm, gds, gmbs;
  if (vgst <= 0.0) {
    gm = 0.0;
    gds = 0.0;
    gmbs = 0.0;
  } else if (vds <= vgst) {
    double k = model.kp;
    gm = k * vds * (1.0 + model.lambda * vds);
    gds = k * (vgst - vds) * (1.0 + model.lambda * vds) + k * (vgst * vds - 0.5 * vds * vds) * model.lambda;
    gmbs = -model.gamma * 0.5 / std::sqrt(std::max(1e-12, model.phi - vbs)) * gm;
  } else {
    double k = model.kp;
    gm = k * vgst * (1.0 + model.lambda * vds);
    gds = 0.5 * k * vgst * vgst * model.lambda;
    gmbs = -model.gamma * 0.5 / std::sqrt(std::max(1e-12, model.phi - vbs)) * gm;
  }

  add_element_c(A, nd, nd, complex_t(gds, 0.0));
  add_element_c(A, nd, ns, complex_t(-gds - gm - gmbs, 0.0));
  add_element_c(A, nd, ng, complex_t(gm, 0.0));
  add_element_c(A, nd, nb, complex_t(gmbs, 0.0));
  add_element_c(A, ns, nd, complex_t(-gds, 0.0));
  add_element_c(A, ns, ns, complex_t(gds + gm + gmbs, 0.0));
  add_element_c(A, ns, ng, complex_t(-gm, 0.0));
  add_element_c(A, ns, nb, complex_t(-gmbs, 0.0));

  double cgs = model.cgso;
  double cgd = model.cgdo;
  complex_t Ygs(0.0, omega * cgs);
  complex_t Ygd(0.0, omega * cgd);
  add_element_c(A, ng, ns, Ygs);
  add_element_c(A, ns, ng, -Ygs);
  add_element_c(A, ng, ng, Ygs + Ygd);
  add_element_c(A, ng, nd, -Ygd);
  add_element_c(A, nd, ng, -Ygd);
  add_element_c(A, nd, nd, Ygd);
}

void stamp_bjt_dc(const Component& comp, const std::map<std::string, int>& node_idx,
                  const BjtModel& model, Eigen::SparseMatrix<double>& A,
                  Eigen::VectorXd& b, const Eigen::VectorXd& x,
                  double vt) {
  int nc = node_idx.at(comp.node_pins[0]);
  int nb = node_idx.at(comp.node_pins[1]);
  int ne = node_idx.at(comp.node_pins[2]);

  double vbe = x(nb) - x(ne);
  double vbc = x(nb) - x(nc);

  double is = model.is;
  double bf = model.bf;
  double nfvt = model.nf * vt;
  double nrvt = model.nr * vt;

  double ic, ib, ie;
  double gpi, gmu, g0, gm_trans, go;

  double vbe_c = std::min(vbe, 100.0 * nfvt);
  double vbc_c = std::min(vbc, 100.0 * nrvt);

  double vbe_min = std::max(vbe_c, -10.0 * nfvt);
  double vbc_min = std::max(vbc_c, -10.0 * nrvt);

  double ifwd = is * (std::exp(vbe_min / nfvt) - 1.0);
  double irev = is * (std::exp(vbc_min / nrvt) - 1.0);

  gm_trans = ifwd / nfvt;

  if (model.vaf > 0.0) {
    go = ifwd / model.vaf;
  } else {
    go = 0.0;
  }

  gpi = gm_trans / bf;
  gmu = 0.0;
  g0 = go;

  ic = ifwd - irev / model.br + g0 * (x(nc) - x(ne));
  ie = -ifwd - ifwd / bf;
  ib = -ic - ie;

  double ic_eq = ic - (gm_trans - gmu) * vbe - gmu * vbc - g0 * (x(nc) - x(ne));
  double ib_eq = ib - gpi * vbe + gmu * vbc;
  double ie_eq = -ic_eq - ib_eq;

  double gbe = gpi + gm_trans;

  add_element(A, nc, nc, g0 + gmu);
  add_element(A, nc, nb, gm_trans - gmu);
  add_element(A, nc, ne, -g0 - gm_trans);
  add_element(A, nb, nc, -gmu);
  add_element(A, nb, nb, gbe + gmu);
  add_element(A, nb, ne, -gbe);
  add_element(A, ne, nc, -g0);
  add_element(A, ne, nb, -gbe - gm_trans);
  add_element(A, ne, ne, g0 + gbe + gm_trans);

  b(nc) -= ic_eq;
  b(nb) -= ib_eq;
  b(ne) -= ie_eq;
}

void stamp_bjt_ac(const Component& comp, const std::map<std::string, int>& node_idx,
                  const BjtModel& model, Eigen::SparseMatrix<complex_t>& A,
                  Eigen::VectorXcd& b, const Eigen::VectorXd& op,
                  double vt, double omega) {
  int nc = node_idx.at(comp.node_pins[0]);
  int nb = node_idx.at(comp.node_pins[1]);
  int ne = node_idx.at(comp.node_pins[2]);

  double vbe = op(nb) - op(ne);
  double nfvt = model.nf * vt;

  double vbe_c = std::min(vbe, 100.0 * nfvt);
  double vbe_min = std::max(vbe_c, -10.0 * nfvt);

  double ifwd = model.is * (std::exp(vbe_min / nfvt) - 1.0);
  double gm_trans = ifwd / nfvt;

  double gpi = gm_trans / model.bf;
  double go = (model.vaf > 0.0) ? ifwd / model.vaf : 0.0;

  double gbe = gpi + gm_trans;

  complex_t Ypi(gpi, omega * model.cje);
  complex_t Ymu(0.0, omega * model.cjc);
  complex_t Ygo(go, 0.0);

  add_element_c(A, nc, nc, Ygo + Ymu);
  add_element_c(A, nc, nb, complex_t(gm_trans, 0.0) - Ymu);
  add_element_c(A, nc, ne, -Ygo - complex_t(gm_trans, 0.0));
  add_element_c(A, nb, nc, -Ymu);
  add_element_c(A, nb, nb, Ypi + Ymu);
  add_element_c(A, nb, ne, -Ypi);
  add_element_c(A, ne, nc, -Ygo);
  add_element_c(A, ne, nb, -Ypi - complex_t(gm_trans, 0.0));
  add_element_c(A, ne, ne, Ygo + Ypi + complex_t(gm_trans, 0.0));
}

}  // namespace ahkab
