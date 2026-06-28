#include "analysis/ac_analysis.hpp"

#include <cmath>

#include "analysis/op_analysis.hpp"
#include "mna/linear_solver.hpp"

namespace ahkab {

AcResult run_ac(const Circuit& circuit, const AcConfig& config) {
  OpResult op = run_op(circuit);

  int N_eff = circuit.node_count_effective();
  int M = circuit.num_vsrcs();
  int K = circuit.num_inductors();
  int dim = N_eff + M + K;

  AcResult result;
  for (const auto& name : circuit.node_names()) {
    result.variable_names.push_back("V(" + name + ")");
  }
  for (const auto& vsrc : circuit.vsrcs()) {
    result.variable_names.push_back("I(" + vsrc.id + ")");
  }

  std::vector<double> freqs;
  if (config.sweep_type == SweepType::LINEAR) {
    double df = (config.fstop - config.fstart) / (config.npts - 1);
    for (int i = 0; i < config.npts; ++i) {
      freqs.push_back(config.fstart + i * df);
    }
  } else {
    double log_start = std::log10(config.fstart);
    double log_stop = std::log10(config.fstop);
    double dlog = (log_stop - log_start) / (config.npts - 1);
    for (int i = 0; i < config.npts; ++i) {
      freqs.push_back(std::pow(10.0, log_start + i * dlog));
    }
  }

  std::map<std::string, int> vsrc_row, ind_row;
  int vr = N_eff;
  for (const auto& vsrc : circuit.vsrcs()) vsrc_row[vsrc.id] = vr++;
  int ir = N_eff + M;
  for (const auto& ind : circuit.inductors()) ind_row[ind.id] = ir++;

  auto ni = [&](const std::string& n) { return circuit.mna_node_index(n); };
  auto add_c = [&](Eigen::SparseMatrix<complex_t>& A, int r, int c, complex_t v) {
    if (r >= 0 && c >= 0) A.coeffRef(r, c) += v;
  };
  auto add_bc = [&](Eigen::VectorXcd& b, int r, complex_t v) {
    if (r >= 0) b(r) += v;
  };

  for (double freq : freqs) {
    double omega = 2.0 * M_PI * freq;
    Eigen::SparseMatrix<complex_t> A_mat(dim, dim);
    Eigen::VectorXcd b_vec(dim);
    b_vec.setZero();

    for (const auto& comp : circuit.components()) {
      switch (comp.type) {
        case ComponentType::RESISTOR: {
          int n1 = ni(comp.node_pins[0]);
          int n2 = ni(comp.node_pins[1]);
          complex_t gz(1.0 / comp.params.at("R"), 0.0);
          add_c(A_mat, n1, n1, gz);
          add_c(A_mat, n1, n2, -gz);
          add_c(A_mat, n2, n1, -gz);
          add_c(A_mat, n2, n2, gz);
          break;
        }
        case ComponentType::CAPACITOR: {
          int n1 = ni(comp.node_pins[0]);
          int n2 = ni(comp.node_pins[1]);
          complex_t Y(0.0, omega * comp.params.at("C"));
          add_c(A_mat, n1, n1, Y);
          add_c(A_mat, n1, n2, -Y);
          add_c(A_mat, n2, n1, -Y);
          add_c(A_mat, n2, n2, Y);
          break;
        }
        case ComponentType::INDUCTOR: {
          int n1 = ni(comp.node_pins[0]);
          int n2 = ni(comp.node_pins[1]);
          int ir_idx = ind_row.at(comp.id);
          complex_t Z(0.0, omega * comp.params.at("L"));
          add_c(A_mat, n1, ir_idx, 1.0);
          add_c(A_mat, n2, ir_idx, -1.0);
          add_c(A_mat, ir_idx, n1, 1.0);
          add_c(A_mat, ir_idx, n2, -1.0);
          add_c(A_mat, ir_idx, ir_idx, -Z);
          break;
        }
        case ComponentType::VSOURCE: {
          int n_pos = ni(comp.node_pins[0]);
          int n_neg = ni(comp.node_pins[1]);
          int vr_idx = vsrc_row.at(comp.id);
          double phase_rad = comp.ac_phase * M_PI / 180.0;
          complex_t val(comp.ac_magnitude * std::cos(phase_rad),
                        comp.ac_magnitude * std::sin(phase_rad));
          add_c(A_mat, n_pos, vr_idx, 1.0);
          add_c(A_mat, n_neg, vr_idx, -1.0);
          add_c(A_mat, vr_idx, n_pos, 1.0);
          add_c(A_mat, vr_idx, n_neg, -1.0);
          b_vec(vr_idx) = val;
          break;
        }
        case ComponentType::ISOURCE: {
          int n_pos = ni(comp.node_pins[0]);
          int n_neg = ni(comp.node_pins[1]);
          double phase_rad = comp.ac_phase * M_PI / 180.0;
          complex_t val(comp.ac_magnitude * std::cos(phase_rad),
                        comp.ac_magnitude * std::sin(phase_rad));
          add_bc(b_vec, n_pos, -val);
          add_bc(b_vec, n_neg, val);
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
          auto xn = [&](int idx) { return idx >= 0 ? op.meta.x(idx) : 0.0; };
          double vgs = xn(ng) - xn(ns);
          double vds = xn(nd) - xn(ns);
          double vbs = xn(nb) - xn(ns);
          double vth = mdl->vto + mdl->gamma *
            (std::sqrt(std::max(0.0, mdl->phi - vbs)) - std::sqrt(mdl->phi));
          double vgst = vgs - vth;
          double gm = 0.0, gds = 0.0, gmbs = 0.0;
          if (vgst > 0.0) {
            if (vds <= vgst) {
              double k = mdl->kp;
              gm = k * vds * (1.0 + mdl->lambda * vds);
              gds = k * (vgst - vds) * (1.0 + mdl->lambda * vds) +
                    k * (vgst * vds - 0.5 * vds * vds) * mdl->lambda;
            } else {
              double k = mdl->kp;
              gm = k * vgst * (1.0 + mdl->lambda * vds);
              gds = 0.5 * k * vgst * vgst * mdl->lambda;
            }
            gmbs = -mdl->gamma * 0.5 / std::sqrt(std::max(1e-12, mdl->phi - vbs)) * gm;
          }
          add_c(A_mat, nd, nd, complex_t(gds, 0.0));
          add_c(A_mat, nd, ns, complex_t(-gds - gm - gmbs, 0.0));
          add_c(A_mat, nd, ng, complex_t(gm, 0.0));
          add_c(A_mat, nd, nb, complex_t(gmbs, 0.0));
          add_c(A_mat, ns, nd, complex_t(-gds, 0.0));
          add_c(A_mat, ns, ns, complex_t(gds + gm + gmbs, 0.0));
          add_c(A_mat, ns, ng, complex_t(-gm, 0.0));
          add_c(A_mat, ns, nb, complex_t(-gmbs, 0.0));
          break;
        }
        case ComponentType::BJT_NPN:
        case ComponentType::BJT_PNP: {
          auto* mdl = circuit.find_bjt_model(comp.model_name);
          if (!mdl) break;
          int nc = ni(comp.node_pins[0]);
          int nb = ni(comp.node_pins[1]);
          int ne = ni(comp.node_pins[2]);
          auto xn = [&](int idx) { return idx >= 0 ? op.meta.x(idx) : 0.0; };
          double vbe = xn(nb) - xn(ne);
          double nfvt = mdl->nf * Vt_thermal();
          double vbe_c = std::min(vbe, 100.0 * nfvt);
          double vbe_min = std::max(vbe_c, -10.0 * nfvt);
          double ifwd = mdl->is * (std::exp(vbe_min / nfvt) - 1.0);
          double gm_trans = ifwd / nfvt;
          double gpi = gm_trans / mdl->bf;
          double go = (mdl->vaf > 0.0) ? ifwd / mdl->vaf : 0.0;
          add_c(A_mat, nc, nc, complex_t(go, 0.0));
          add_c(A_mat, nc, nb, complex_t(gm_trans, 0.0));
          add_c(A_mat, nc, ne, complex_t(-go - gm_trans, 0.0));
          add_c(A_mat, nb, nb, complex_t(gpi, 0.0));
          add_c(A_mat, nb, ne, complex_t(-gpi, 0.0));
          add_c(A_mat, ne, nc, complex_t(-go, 0.0));
          add_c(A_mat, ne, nb, complex_t(-gpi - gm_trans, 0.0));
          add_c(A_mat, ne, ne, complex_t(go + gpi + gm_trans, 0.0));
          break;
        }
        default:
          break;
      }
    }

    A_mat.makeCompressed();
    Eigen::VectorXcd sol = solve_linear_complex(A_mat, b_vec);

    result.frequencies.push_back(freq);
    std::vector<complex_t> row;
    for (const auto& name : circuit.node_names()) {
      int idx = ni(name);
      row.push_back(idx >= 0 ? sol(idx) : complex_t(0.0, 0.0));
    }
    for (int i = N_eff; i < dim; ++i) {
      row.push_back(sol(i));
    }
    result.values.push_back(row);
  }

  return result;
}

}  // namespace ahkab
