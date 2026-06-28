#include "mna/nr_solver.hpp"

#include <cmath>

#include "core/errors.hpp"
#include "mna/linear_solver.hpp"
#include "mna/mna_builder.hpp"

namespace ahkab {

NrResult solve_dc(const Circuit& circuit, const NrConfig& config) {
  int dim = circuit.mna_dimension();
  NrResult result;
  result.x = Eigen::VectorXd::Zero(dim);

  auto ni = [&](const std::string& n) { return circuit.mna_node_index(n); };
  auto add = [&](Eigen::SparseMatrix<double>& A, int r, int c, double v) {
    if (r >= 0 && c >= 0) A.coeffRef(r, c) += v;
  };
  auto add_b = [&](Eigen::VectorXd& b, int r, double v) {
    if (r >= 0) b(r) += v;
  };

  auto build_sys = [&](const Eigen::VectorXd& x) -> MnaSystem {
    MnaSystem sys = build_mna_dc(circuit);
    for (const auto& comp : circuit.components()) {
      switch (comp.type) {
        case ComponentType::DIODE: {
          auto* model = circuit.find_diode_model(comp.model_name);
          if (model) {
            int na = ni(comp.node_pins[0]);
            int nc = ni(comp.node_pins[1]);
            double vd = 0.0;
            if (na >= 0 && nc >= 0) vd = x(na) - x(nc);
            else if (na >= 0) vd = x(na);
            else if (nc >= 0) vd = -x(nc);
            double nvt = model->n * Vt_thermal();
            double vcrit = nvt * std::log(nvt / (std::sqrt(2.0) * model->is));
            double id, gd;
            if (vd >= vcrit + 10.0 * nvt) {
              double evd = std::exp(vcrit / nvt);
              id = model->is * evd * (1.0 + (vd - vcrit) / nvt);
              gd = model->is * evd / nvt;
            } else if (vd < -5.0 * nvt) {
              id = -model->is;
              gd = 0.0;
            } else {
              double evd = std::exp(vd / nvt);
              id = model->is * (evd - 1.0);
              gd = model->is * evd / nvt;
            }
            double ieq = id - gd * vd;
            add(sys.A, na, na, gd);
            add(sys.A, na, nc, -gd);
            add(sys.A, nc, na, -gd);
            add(sys.A, nc, nc, gd);
            add_b(sys.b, na, -ieq);
            add_b(sys.b, nc, ieq);
          }
          break;
        }
        case ComponentType::MOSFET_N:
        case ComponentType::MOSFET_P: {
          auto* model = circuit.find_mosfet_model(comp.model_name);
          if (model) {
            int nd = ni(comp.node_pins[0]);
            int ng = ni(comp.node_pins[1]);
            int ns = ni(comp.node_pins[2]);
            int nb = ni(comp.node_pins[3]);
            auto xn = [&](int idx) { return idx >= 0 ? x(idx) : 0.0; };
            double vgs = xn(ng) - xn(ns);
            double vds = xn(nd) - xn(ns);
            double vbs = xn(nb) - xn(ns);
            double vth = model->vto + model->gamma *
              (std::sqrt(std::max(0.0, model->phi - vbs)) - std::sqrt(model->phi));
            double vgst = vgs - vth;
            double ids, gm, gds, gmbs;
            if (vgst <= 0.0) {
              ids = 0.0; gm = 0.0; gds = 0.0; gmbs = 0.0;
            } else if (vds <= vgst) {
              double k = model->kp;
              ids = k * (vgst * vds - 0.5 * vds * vds) * (1.0 + model->lambda * vds);
              gm = k * vds * (1.0 + model->lambda * vds);
              gds = k * (vgst - vds) * (1.0 + model->lambda * vds) +
                    k * (vgst * vds - 0.5 * vds * vds) * model->lambda;
              gmbs = -model->gamma * 0.5 / std::sqrt(std::max(1e-12, model->phi - vbs)) * gm;
            } else {
              double k = model->kp;
              ids = 0.5 * k * vgst * vgst * (1.0 + model->lambda * vds);
              gm = k * vgst * (1.0 + model->lambda * vds);
              gds = 0.5 * k * vgst * vgst * model->lambda;
              gmbs = -model->gamma * 0.5 / std::sqrt(std::max(1e-12, model->phi - vbs)) * gm;
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
          }
          break;
        }
        case ComponentType::BJT_NPN:
        case ComponentType::BJT_PNP: {
          auto* model = circuit.find_bjt_model(comp.model_name);
          if (model) {
            int nc = ni(comp.node_pins[0]);
            int nb = ni(comp.node_pins[1]);
            int ne = ni(comp.node_pins[2]);
            auto xn = [&](int idx) { return idx >= 0 ? x(idx) : 0.0; };
            double vt = Vt_thermal();
            double vbe = xn(nb) - xn(ne);
            double nfvt = model->nf * vt;
            double vbe_c = std::min(vbe, 100.0 * nfvt);
            double vbe_min = std::max(vbe_c, -10.0 * nfvt);
            double ifwd = model->is * (std::exp(vbe_min / nfvt) - 1.0);
            double gm_trans = ifwd / nfvt;
            double gpi = gm_trans / model->bf;
            double go = (model->vaf > 0.0) ? ifwd / model->vaf : 0.0;
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
          }
          break;
        }
        default:
          break;
      }
    }
    return sys;
  };

  double current_gmin = config.gmin;
  for (int gs = config.gmin_steps; gs >= 0; --gs) {
    if (gs == 0) {
      current_gmin = 0.0;
    } else {
      current_gmin = config.gmin * std::pow(10.0, gs - config.gmin_steps);
    }

    for (int iter = 0; iter < config.max_iter; ++iter) {
      MnaSystem sys = build_sys(result.x);

      int N_eff = circuit.node_count_effective();
      for (int i = 0; i < N_eff; ++i) {
        for (int j = 0; j < N_eff; ++j) {
          if (i != j) {
            sys.A.coeffRef(i, j) += current_gmin * (-1.0 / N_eff);
          }
        }
        sys.A.coeffRef(i, i) += current_gmin;
      }

      sys.A.makeCompressed();
      Eigen::VectorXd r = sys.b - sys.A * result.x;
      Eigen::VectorXd dx = solve_linear(sys.A, r);

      result.x += dx;
      result.iterations++;

      double rnorm = dx.lpNorm<Eigen::Infinity>();
      double xnorm = result.x.lpNorm<Eigen::Infinity>();
      result.residual = rnorm;

      if (rnorm < config.abstol * dim + config.reltol * xnorm) {
        result.converged = true;
        return result;
      }
    }
  }

  throw ConvergenceError("DC solve failed to converge",
                          result.residual, result.iterations);
}

}  // namespace ahkab
