#include "core/circuit.hpp"

#include <stdexcept>

namespace ahkab {

int Circuit::register_or_get_node(const std::string& name) {
  auto it = m_node_index.find(name);
  if (it != m_node_index.end()) {
    return it->second;
  }
  int idx = static_cast<int>(m_node_names.size());
  m_node_names.push_back(name);
  m_node_index[name] = idx;
  return idx;
}

int Circuit::node_count() const {
  return static_cast<int>(m_node_names.size());
}

int Circuit::node_count_effective() const {
  if (m_ground.empty()) return node_count();
  return node_count() - 1;
}

int Circuit::num_vsrcs() const {
  return static_cast<int>(m_vsrcs.size());
}

int Circuit::num_inductors() const {
  return static_cast<int>(m_inductors.size());
}

int Circuit::mna_dimension() const {
  return node_count_effective() + num_vsrcs() + num_inductors();
}

int Circuit::mna_node_index(const std::string& name) const {
  auto it = m_node_index.find(name);
  if (it == m_node_index.end()) {
    throw std::runtime_error("Node not found: " + name);
  }
  int raw = it->second;
  int gnd = m_node_index.at(m_ground);
  if (raw == gnd) return -1;
  if (raw > gnd) return raw - 1;
  return raw;
}

void Circuit::set_ground(const std::string& node_name) {
  m_ground = node_name;
  register_or_get_node(node_name);
}

int Circuit::ground_index() const {
  auto it = m_node_index.find(m_ground);
  if (it == m_node_index.end()) {
    throw std::runtime_error("Ground node not registered");
  }
  return it->second;
}

void Circuit::add_component(const Component& comp) {
  m_components.push_back(comp);
  if (comp.type == ComponentType::VSOURCE ||
      comp.type == ComponentType::VCVS ||
      comp.type == ComponentType::CCVS) {
    m_vsrcs.push_back(comp);
  }
  if (comp.type == ComponentType::INDUCTOR) {
    m_inductors.push_back(comp);
  }
}

void Circuit::add_resistor(const std::string& id, const std::string& n1,
                           const std::string& n2, double value) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  Component c;
  c.id = id;
  c.type = ComponentType::RESISTOR;
  c.node_pins = {n1, n2};
  c.params["R"] = value;
  add_component(c);
}

void Circuit::add_capacitor(const std::string& id, const std::string& n1,
                            const std::string& n2, double value, double ic) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  Component c;
  c.id = id;
  c.type = ComponentType::CAPACITOR;
  c.node_pins = {n1, n2};
  c.params["C"] = value;
  if (ic != 0.0) {
    m_ic["V(" + n1 + "," + n2 + ")"] = ic;
  }
  add_component(c);
}

void Circuit::add_inductor(const std::string& id, const std::string& n1,
                           const std::string& n2, double value, double ic) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  Component c;
  c.id = id;
  c.type = ComponentType::INDUCTOR;
  c.node_pins = {n1, n2};
  c.params["L"] = value;
  if (ic != 0.0) {
    c.dc_value = ic;
  }
  add_component(c);
}

void Circuit::add_vsource(const std::string& id, const std::string& n_pos,
                          const std::string& n_neg, double dc_value,
                          double ac_mag, double ac_phase) {
  register_or_get_node(n_pos);
  register_or_get_node(n_neg);
  Component c;
  c.id = id;
  c.type = ComponentType::VSOURCE;
  c.node_pins = {n_pos, n_neg};
  c.dc_value = dc_value;
  c.ac_magnitude = ac_mag;
  c.ac_phase = ac_phase;
  add_component(c);
}

void Circuit::add_isource(const std::string& id, const std::string& n_pos,
                          const std::string& n_neg, double dc_value,
                          double ac_mag, double ac_phase) {
  register_or_get_node(n_pos);
  register_or_get_node(n_neg);
  Component c;
  c.id = id;
  c.type = ComponentType::ISOURCE;
  c.node_pins = {n_pos, n_neg};
  c.dc_value = dc_value;
  c.ac_magnitude = ac_mag;
  c.ac_phase = ac_phase;
  add_component(c);
}

void Circuit::add_vcvs(const std::string& id, const std::string& n1,
                       const std::string& n2, const std::string& ctrl_n1,
                       const std::string& ctrl_n2, double gain) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  register_or_get_node(ctrl_n1);
  register_or_get_node(ctrl_n2);
  Component c;
  c.id = id;
  c.type = ComponentType::VCVS;
  c.node_pins = {n1, n2, ctrl_n1, ctrl_n2};
  c.params["gain"] = gain;
  add_component(c);
}

void Circuit::add_vccs(const std::string& id, const std::string& n1,
                       const std::string& n2, const std::string& ctrl_n1,
                       const std::string& ctrl_n2, double gain) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  register_or_get_node(ctrl_n1);
  register_or_get_node(ctrl_n2);
  Component c;
  c.id = id;
  c.type = ComponentType::VCCS;
  c.node_pins = {n1, n2, ctrl_n1, ctrl_n2};
  c.params["gain"] = gain;
  add_component(c);
}

void Circuit::add_ccvs(const std::string& id, const std::string& n1,
                       const std::string& n2, const std::string& vbranch,
                       double gain) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  Component c;
  c.id = id;
  c.type = ComponentType::CCVS;
  c.node_pins = {n1, n2, vbranch};
  c.params["gain"] = gain;
  add_component(c);
}

void Circuit::add_cccs(const std::string& id, const std::string& n1,
                       const std::string& n2, const std::string& vbranch,
                       double gain) {
  register_or_get_node(n1);
  register_or_get_node(n2);
  Component c;
  c.id = id;
  c.type = ComponentType::CCCS;
  c.node_pins = {n1, n2, vbranch};
  c.params["gain"] = gain;
  add_component(c);
}

void Circuit::add_diode(const std::string& id, const std::string& anode,
                        const std::string& cathode,
                        const std::string& model_name) {
  register_or_get_node(anode);
  register_or_get_node(cathode);
  Component c;
  c.id = id;
  c.type = ComponentType::DIODE;
  c.node_pins = {anode, cathode};
  c.model_name = model_name;
  add_component(c);
}

void Circuit::add_mosfet(const std::string& id, const std::string& drain,
                         const std::string& gate, const std::string& source,
                         const std::string& bulk, const std::string& model_name,
                         ComponentType type) {
  register_or_get_node(drain);
  register_or_get_node(gate);
  register_or_get_node(source);
  register_or_get_node(bulk);
  Component c;
  c.id = id;
  c.type = type;
  c.node_pins = {drain, gate, source, bulk};
  c.model_name = model_name;
  add_component(c);
}

void Circuit::add_bjt(const std::string& id, const std::string& collector,
                      const std::string& base, const std::string& emitter,
                      const std::string& model_name, ComponentType type) {
  register_or_get_node(collector);
  register_or_get_node(base);
  register_or_get_node(emitter);
  Component c;
  c.id = id;
  c.type = type;
  c.node_pins = {collector, base, emitter};
  c.model_name = model_name;
  add_component(c);
}

void Circuit::add_diode_model(const std::string& name,
                               const DiodeModel& model) {
  DiodeModel m = model;
  m.name = name;
  m_diode_models[name] = m;
}

void Circuit::add_mosfet_model(const std::string& name,
                                const MosfetModel& model) {
  MosfetModel m = model;
  m.name = name;
  m_mosfet_models[name] = m;
}

void Circuit::add_bjt_model(const std::string& name, const BjtModel& model) {
  BjtModel m = model;
  m.name = name;
  m_bjt_models[name] = m;
}

const DiodeModel* Circuit::find_diode_model(const std::string& name) const {
  auto it = m_diode_models.find(name);
  return (it != m_diode_models.end()) ? &it->second : nullptr;
}

const MosfetModel* Circuit::find_mosfet_model(const std::string& name) const {
  auto it = m_mosfet_models.find(name);
  return (it != m_mosfet_models.end()) ? &it->second : nullptr;
}

const BjtModel* Circuit::find_bjt_model(const std::string& name) const {
  auto it = m_bjt_models.find(name);
  return (it != m_bjt_models.end()) ? &it->second : nullptr;
}

}  // namespace ahkab
