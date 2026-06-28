#include "mna/mna_builder.hpp"

namespace ahkab {

MnaSystem build_mna_dc(const Circuit& circuit) {
  MnaSystem sys;
  int N_eff = circuit.node_count_effective();
  int M = circuit.num_vsrcs();
  int K = circuit.num_inductors();
  sys.dimension = N_eff + M + K;
  sys.A.resize(sys.dimension, sys.dimension);
  sys.b.resize(sys.dimension);
  sys.b.setZero();

  int vsrc_counter = N_eff;
  for (const auto& vsrc : circuit.vsrcs()) {
    sys.vsrc_row[vsrc.id] = vsrc_counter++;
  }

  int ind_counter = N_eff + M;
  for (const auto& ind : circuit.inductors()) {
    sys.ind_row[ind.id] = ind_counter++;
  }

  auto ni = [&](const std::string& n) { return circuit.mna_node_index(n); };
  auto vr = [&](const std::string& id) { return sys.vsrc_row.at(id); };
  auto ir = [&](const std::string& id) { return sys.ind_row.at(id); };

  auto add = [&](int r, int c, double v) {
    if (r >= 0 && c >= 0) sys.A.coeffRef(r, c) += v;
  };
  auto add_b = [&](int r, double v) {
    if (r >= 0) sys.b(r) += v;
  };

  for (const auto& comp : circuit.components()) {
    switch (comp.type) {
      case ComponentType::RESISTOR: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        double g = 1.0 / comp.params.at("R");
        add(n1, n1, g);
        add(n1, n2, -g);
        add(n2, n1, -g);
        add(n2, n2, g);
        break;
      }
      case ComponentType::CAPACITOR: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        double gm = 1e-12;
        add(n1, n1, gm);
        add(n1, n2, -gm);
        add(n2, n1, -gm);
        add(n2, n2, gm);
        break;
      }
      case ComponentType::INDUCTOR: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        int ind_r = ir(comp.id);
        double gm = 1e-12;
        add(n1, n1, gm);
        add(n1, n2, -gm);
        add(n2, n1, -gm);
        add(n2, n2, gm);
        add(n1, ind_r, 1.0);
        add(n2, ind_r, -1.0);
        add(ind_r, n1, 1.0);
        add(ind_r, n2, -1.0);
        break;
      }
      case ComponentType::VSOURCE: {
        int n_pos = ni(comp.node_pins[0]);
        int n_neg = ni(comp.node_pins[1]);
        int vsrc_r = vr(comp.id);
        add(n_pos, vsrc_r, 1.0);
        add(n_neg, vsrc_r, -1.0);
        add(vsrc_r, n_pos, 1.0);
        add(vsrc_r, n_neg, -1.0);
        sys.b(vsrc_r) = comp.dc_value;
        break;
      }
      case ComponentType::ISOURCE: {
        int n_pos = ni(comp.node_pins[0]);
        int n_neg = ni(comp.node_pins[1]);
        add_b(n_pos, -comp.dc_value);
        add_b(n_neg, comp.dc_value);
        break;
      }
      case ComponentType::VCVS: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        int cn1 = ni(comp.node_pins[2]);
        int cn2 = ni(comp.node_pins[3]);
        int vsrc_r = vr(comp.id);
        double gain = comp.params.at("gain");
        add(n1, vsrc_r, 1.0);
        add(n2, vsrc_r, -1.0);
        add(vsrc_r, n1, 1.0);
        add(vsrc_r, n2, -1.0);
        add(vsrc_r, cn1, -gain);
        add(vsrc_r, cn2, gain);
        break;
      }
      case ComponentType::VCCS: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        int cn1 = ni(comp.node_pins[2]);
        int cn2 = ni(comp.node_pins[3]);
        double gm = comp.params.at("gain");
        add(n1, cn1, gm);
        add(n1, cn2, -gm);
        add(n2, cn1, -gm);
        add(n2, cn2, gm);
        break;
      }
      case ComponentType::CCVS: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        int vbr = vr(comp.node_pins[2]);
        int vsrc_r = vr(comp.id);
        double gain = comp.params.at("gain");
        add(n1, vsrc_r, 1.0);
        add(n2, vsrc_r, -1.0);
        add(vsrc_r, n1, 1.0);
        add(vsrc_r, n2, -1.0);
        add(vsrc_r, vbr, -gain);
        break;
      }
      case ComponentType::CCCS: {
        int n1 = ni(comp.node_pins[0]);
        int n2 = ni(comp.node_pins[1]);
        int vbr = vr(comp.node_pins[2]);
        double f = comp.params.at("gain");
        add(n1, vbr, f);
        add(n2, vbr, -f);
        break;
      }
      default:
        break;
    }
  }

  return sys;
}

}  // namespace ahkab
