#include "analysis/transient_analysis.hpp"

#include <cmath>

#include "analysis/op_analysis.hpp"
#include "core/errors.hpp"
#include "mna/linear_solver.hpp"
#include "mna/mna_builder.hpp"

namespace ahkab {

TranResult run_tran(const Circuit& circuit, const TranConfig& config) {
  int N_eff = circuit.node_count_effective();
  int M = circuit.num_vsrcs();
  int K = circuit.num_inductors();
  int dim = N_eff + M + K;

  OpResult op = run_op(circuit, config.nr);
  Eigen::VectorXd x = op.meta.x;
  Eigen::VectorXd x_prev = x;
  double t = config.tstart;

  TranResult result;
  result.time.push_back(t);

  for (const auto& name : circuit.node_names()) {
    result.variable_names.push_back("V(" + name + ")");
  }
  for (const auto& vsrc : circuit.vsrcs()) {
    result.variable_names.push_back("I(" + vsrc.id + ")");
  }

  // Map MNA solution back to all nodes including ground
  auto map_to_nodes = [&](const Eigen::VectorXd& sol) -> std::vector<double> {
    std::vector<double> row;
    for (const auto& name : circuit.node_names()) {
      int idx = circuit.mna_node_index(name);
      row.push_back(idx >= 0 ? sol(idx) : 0.0);
    }
    for (int i = N_eff; i < dim; ++i) {
      row.push_back(sol(i));
    }
    return row;
  };

  result.values.push_back(map_to_nodes(x));

  auto ni = [&](const std::string& n) { return circuit.mna_node_index(n); };
  auto add = [&](Eigen::SparseMatrix<double>& A, int r, int c, double v) {
    if (r >= 0 && c >= 0) A.coeffRef(r, c) += v;
  };
  auto add_b = [&](Eigen::VectorXd& b, int r, double v) {
    if (r >= 0) b(r) += v;
  };

  std::map<std::string, double> cap_prev_current;
  for (const auto& comp : circuit.components()) {
    if (comp.type == ComponentType::CAPACITOR)
      cap_prev_current[comp.id] = 0.0;
  }

  double dt = config.tstep;

  while (t < config.tstop - dt * 0.5) {
    double dt_eff = std::min(dt, config.tmax);
    bool step_accepted = false;
    Eigen::VectorXd x_new;

    while (!step_accepted && dt_eff > config.tstep * 1e-12) {
      MnaSystem sys;
      sys.dimension = dim;
      sys.A.resize(dim, dim);
      sys.b.resize(dim);
      sys.b.setZero();

      int vc = N_eff;
      for (const auto& vsrc : circuit.vsrcs()) sys.vsrc_row[vsrc.id] = vc++;
      int ic = N_eff + M;
      for (const auto& ind : circuit.inductors()) sys.ind_row[ind.id] = ic++;

      auto vr = [&](const std::string& id) { return sys.vsrc_row.at(id); };
      auto ir = [&](const std::string& id) { return sys.ind_row.at(id); };

      for (const auto& comp : circuit.components()) {
        switch (comp.type) {
          case ComponentType::RESISTOR: {
            int n1 = ni(comp.node_pins[0]);
            int n2 = ni(comp.node_pins[1]);
            double g = 1.0 / comp.params.at("R");
            add(sys.A, n1, n1, g);
            add(sys.A, n1, n2, -g);
            add(sys.A, n2, n1, -g);
            add(sys.A, n2, n2, g);
            break;
          }
           case ComponentType::CAPACITOR: {
            int n1 = ni(comp.node_pins[0]);
            int n2 = ni(comp.node_pins[1]);
            double C = comp.params.at("C");
            double geq;
            if (config.method == IntegrationMethod::IMPLICIT_EULER) {
              geq = C / dt_eff;
            } else if (config.method == IntegrationMethod::TRAPEZOIDAL) {
              geq = 2.0 * C / dt_eff;
            } else {
              geq = 2.0 * C / dt_eff;
            }
            auto xn = [&](int idx) { return idx >= 0 ? x_prev(idx) : 0.0; };
            double vc_prev = xn(n1) - xn(n2);
            double ic_prev = cap_prev_current[comp.id];
            double ieq;
            if (config.method == IntegrationMethod::IMPLICIT_EULER) {
              ieq = geq * vc_prev;
            } else {
              ieq = geq * vc_prev + ic_prev;
            }
            add(sys.A, n1, n1, geq);
            add(sys.A, n1, n2, -geq);
            add(sys.A, n2, n1, -geq);
            add(sys.A, n2, n2, geq);
            add_b(sys.b, n1, ieq);
            add_b(sys.b, n2, -ieq);
            break;
          }
          case ComponentType::INDUCTOR: {
            int n1 = ni(comp.node_pins[0]);
            int n2 = ni(comp.node_pins[1]);
            int ind_r = ir(comp.id);
            double L = comp.params.at("L");
            double req;
            if (config.method == IntegrationMethod::IMPLICIT_EULER) {
              req = L / dt_eff;
            } else if (config.method == IntegrationMethod::TRAPEZOIDAL) {
              req = 2.0 * L / dt_eff;
            } else {
              req = 2.0 * L / dt_eff;
            }
            double il_prev = x_prev(ind_r);
            add(sys.A, n1, ind_r, 1.0);
            add(sys.A, n2, ind_r, -1.0);
            add(sys.A, ind_r, n1, 1.0);
            add(sys.A, ind_r, n2, -1.0);
            add(sys.A, ind_r, ind_r, -req);
            auto xn = [&](int idx) { return idx >= 0 ? x_prev(idx) : 0.0; };
            if (config.method == IntegrationMethod::IMPLICIT_EULER) {
              sys.b(ind_r) -= req * il_prev;
            } else {
              double vl_prev = xn(n1) - xn(n2);
              sys.b(ind_r) -= req * il_prev + vl_prev;
            }
            break;
          }
          case ComponentType::VSOURCE: {
            int n_pos = ni(comp.node_pins[0]);
            int n_neg = ni(comp.node_pins[1]);
            int vsrc_r = vr(comp.id);
            add(sys.A, n_pos, vsrc_r, 1.0);
            add(sys.A, n_neg, vsrc_r, -1.0);
            add(sys.A, vsrc_r, n_pos, 1.0);
            add(sys.A, vsrc_r, n_neg, -1.0);
            sys.b(vsrc_r) = comp.dc_value;
            break;
          }
          case ComponentType::ISOURCE: {
            int n_pos = ni(comp.node_pins[0]);
            int n_neg = ni(comp.node_pins[1]);
            add_b(sys.b, n_pos, -comp.dc_value);
            add_b(sys.b, n_neg, comp.dc_value);
            break;
          }
          case ComponentType::DIODE: {
            auto* mdl = circuit.find_diode_model(comp.model_name);
            if (!mdl) break;
            int na = ni(comp.node_pins[0]);
            int nc = ni(comp.node_pins[1]);
            auto xn = [&](int idx) { return idx >= 0 ? x(idx) : 0.0; };
            double vd = xn(na) - xn(nc);
            double nvt = mdl->n * Vt_thermal();
            double vcrit = nvt * std::log(nvt / (std::sqrt(2.0) * mdl->is));
            double id, gd;
            if (vd >= vcrit + 10.0 * nvt) {
              double evd = std::exp(vcrit / nvt);
              id = mdl->is * evd * (1.0 + (vd - vcrit) / nvt);
              gd = mdl->is * evd / nvt;
            } else if (vd < -5.0 * nvt) {
              id = -mdl->is; gd = 0.0;
            } else {
              double evd = std::exp(vd / nvt);
              id = mdl->is * (evd - 1.0);
              gd = mdl->is * evd / nvt;
            }
            double ieq = id - gd * vd;
            add(sys.A, na, na, gd);
            add(sys.A, na, nc, -gd);
            add(sys.A, nc, na, -gd);
            add(sys.A, nc, nc, gd);
            add_b(sys.b, na, -ieq);
            add_b(sys.b, nc, ieq);
            break;
          }
          case ComponentType::MOSFET_N:
          case ComponentType::MOSFET_P: {
            auto* mdl = circuit.find_mosfet_model(comp.model_name);
            if (!mdl) break;
            int nd = ni(comp.node_pins[0]);
            int ng = ni(comp.node_pins[1]);
            int ns = ni(comp.node_pins[2]);
            int nb = ni(comp.node_pins[3]);
            auto xn = [&](int idx) { return idx >= 0 ? x(idx) : 0.0; };
            double vgs = xn(ng) - xn(ns);
            double vds = xn(nd) - xn(ns);
            double vbs = xn(nb) - xn(ns);
            double vth = mdl->vto + mdl->gamma *
              (std::sqrt(std::max(0.0, mdl->phi - vbs)) - std::sqrt(mdl->phi));
            double vgst = vgs - vth;
            double ids, gm, gds, gmbs;
            if (vgst <= 0.0) {
              ids = 0.0; gm = 0.0; gds = 0.0; gmbs = 0.0;
            } else if (vds <= vgst) {
              double k = mdl->kp;
              ids = k * (vgst * vds - 0.5 * vds * vds) * (1.0 + mdl->lambda * vds);
              gm = k * vds * (1.0 + mdl->lambda * vds);
              gds = k * (vgst - vds) * (1.0 + mdl->lambda * vds) +
                    k * (vgst * vds - 0.5 * vds * vds) * mdl->lambda;
              gmbs = -mdl->gamma * 0.5 / std::sqrt(std::max(1e-12, mdl->phi - vbs)) * gm;
            } else {
              double k = mdl->kp;
              ids = 0.5 * k * vgst * vgst * (1.0 + mdl->lambda * vds);
              gm = k * vgst * (1.0 + mdl->lambda * vds);
              gds = 0.5 * k * vgst * vgst * mdl->lambda;
              gmbs = -mdl->gamma * 0.5 / std::sqrt(std::max(1e-12, mdl->phi - vbs)) * gm;
            }
            double ids_eq = ids - gm * vgs - gds * vds - gmbs * vbs;
            add(sys.A, nd, nd, gds);
            add(sys.A, nd, ns, -gds - gm - gmbs);
            add(sys.A, nd, ng, gm);
            add(sys.A, nd, nb, gmbs);
            add(sys.A, ns, nd, -gds);
            add(sys.A, ns, ns, gds + gm + gmbs);
            add(sys.A, ns, ng, -gm);
            add(sys.A, ns, nb, -gmbs);
            add_b(sys.b, nd, -ids_eq);
            add_b(sys.b, ns, ids_eq);
            break;
          }
          case ComponentType::BJT_NPN:
          case ComponentType::BJT_PNP: {
            auto* mdl = circuit.find_bjt_model(comp.model_name);
            if (!mdl) break;
            int nc = ni(comp.node_pins[0]);
            int nb = ni(comp.node_pins[1]);
            int ne = ni(comp.node_pins[2]);
            auto xn = [&](int idx) { return idx >= 0 ? x(idx) : 0.0; };
            double vt = Vt_thermal();
            double vbe = xn(nb) - xn(ne);
            double nfvt = mdl->nf * vt;
            double vbe_c = std::min(vbe, 100.0 * nfvt);
            double vbe_min = std::max(vbe_c, -10.0 * nfvt);
            double ifwd = mdl->is * (std::exp(vbe_min / nfvt) - 1.0);
            double gm_trans = ifwd / nfvt;
            double gpi = gm_trans / mdl->bf;
            double go = (mdl->vaf > 0.0) ? ifwd / mdl->vaf : 0.0;
            double ic_eq = -gm_trans * vbe - go * (xn(nc) - xn(ne));
            double ib_eq = -gpi * vbe;
            double ie_eq = -ic_eq - ib_eq;
            add(sys.A, nc, nc, go);
            add(sys.A, nc, nb, gm_trans);
            add(sys.A, nc, ne, -go - gm_trans);
            add(sys.A, nb, nb, gpi);
            add(sys.A, nb, ne, -gpi);
            add(sys.A, ne, nc, -go);
            add(sys.A, ne, nb, -gpi - gm_trans);
            add(sys.A, ne, ne, go + gpi + gm_trans);
            add_b(sys.b, nc, -ic_eq);
            add_b(sys.b, nb, -ib_eq);
            add_b(sys.b, ne, -ie_eq);
            break;
          }
          default:
            break;
        }
      }

      sys.A.makeCompressed();
      Eigen::VectorXd r = sys.b - sys.A * x;
      Eigen::VectorXd dx = solve_linear(sys.A, r);
      x_new = x + dx;

      double rnorm = dx.lpNorm<Eigen::Infinity>();
      double xnorm = x_new.lpNorm<Eigen::Infinity>();
      if (rnorm < config.nr.abstol * dim + config.nr.reltol * xnorm) {
        step_accepted = true;
      } else {
        dt_eff *= 0.5;
      }
    }

    if (!step_accepted) {
      throw TimeStepError("Time step underflow", t);
    }

    x_prev = x;
    x = x_new;
    t += dt_eff;

    for (const auto& comp : circuit.components()) {
      if (comp.type == ComponentType::CAPACITOR) {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        double C = comp.params.at("C");
        auto xn = [&](int idx) { return idx >= 0 ? x(idx) : 0.0; };
        auto xn_prev = [&](int idx) { return idx >= 0 ? x_prev(idx) : 0.0; };
        double vc_new = xn(n1) - xn(n2);
        double vc_old = xn_prev(n1) - xn_prev(n2);
        double Cdt = C / dt_eff;
        double ic_old = cap_prev_current[comp.id];
        if (config.method == IntegrationMethod::IMPLICIT_EULER) {
          cap_prev_current[comp.id] = Cdt * (vc_new - vc_old);
        } else {
          cap_prev_current[comp.id] = 2.0 * Cdt * (vc_new - vc_old) - ic_old;
        }
      }
    }

    result.time.push_back(t);
    result.values.push_back(map_to_nodes(x));
  }

  return result;
}

}  // namespace ahkab
