#pragma once

#include <map>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace ahkab {

class Circuit {
 public:
  void add_resistor(const std::string& id, const std::string& n1,
                    const std::string& n2, double value);
  void add_capacitor(const std::string& id, const std::string& n1,
                     const std::string& n2, double value, double ic = 0.0);
  void add_inductor(const std::string& id, const std::string& n1,
                    const std::string& n2, double value, double ic = 0.0);
  void add_vsource(const std::string& id, const std::string& n_pos,
                   const std::string& n_neg, double dc_value,
                   double ac_mag = 0.0, double ac_phase = 0.0);
  void add_isource(const std::string& id, const std::string& n_pos,
                   const std::string& n_neg, double dc_value,
                   double ac_mag = 0.0, double ac_phase = 0.0);
  void add_vcvs(const std::string& id, const std::string& n1,
                const std::string& n2, const std::string& ctrl_n1,
                const std::string& ctrl_n2, double gain);
  void add_vccs(const std::string& id, const std::string& n1,
                const std::string& n2, const std::string& ctrl_n1,
                const std::string& ctrl_n2, double gain);
  void add_ccvs(const std::string& id, const std::string& n1,
                const std::string& n2, const std::string& vbranch,
                double gain);
  void add_cccs(const std::string& id, const std::string& n1,
                const std::string& n2, const std::string& vbranch,
                double gain);
  void add_diode(const std::string& id, const std::string& anode,
                 const std::string& cathode, const std::string& model_name);
  void add_mosfet(const std::string& id, const std::string& drain,
                  const std::string& gate, const std::string& source,
                  const std::string& bulk, const std::string& model_name,
                  ComponentType type);
  void add_bjt(const std::string& id, const std::string& collector,
               const std::string& base, const std::string& emitter,
               const std::string& model_name, ComponentType type);

  void add_diode_model(const std::string& name, const DiodeModel& model);
  void add_mosfet_model(const std::string& name, const MosfetModel& model);
  void add_bjt_model(const std::string& name, const BjtModel& model);

  void add_component(const Component& comp);

  int register_or_get_node(const std::string& name);
  int node_count() const;
  int node_count_effective() const;
  int num_vsrcs() const;
  int num_inductors() const;
  int mna_dimension() const;
  int mna_node_index(const std::string& name) const;

  const std::vector<std::string>& node_names() const { return m_node_names; }
  const std::map<std::string, int>& node_index() const { return m_node_index; }
  const std::vector<Component>& components() const { return m_components; }
  std::vector<Component>& components_mut() { return m_components; }
  const std::vector<Component>& vsrcs() const { return m_vsrcs; }
  const std::vector<Component>& inductors() const { return m_inductors; }
  const std::map<std::string, DiodeModel>& diode_models() const {
    return m_diode_models;
  }
  const std::map<std::string, MosfetModel>& mosfet_models() const {
    return m_mosfet_models;
  }
  const std::map<std::string, BjtModel>& bjt_models() const {
    return m_bjt_models;
  }

  const DiodeModel* find_diode_model(const std::string& name) const;
  const MosfetModel* find_mosfet_model(const std::string& name) const;
  const BjtModel* find_bjt_model(const std::string& name) const;

  void set_ground(const std::string& node_name);
  std::string ground_node() const { return m_ground; }
  int ground_index() const;

  bool has_ground() const { return !m_ground.empty(); }
  const std::map<std::string, double>& initial_conditions() const {
    return m_ic;
  }

  void set_title(const std::string& t) { m_title = t; }
  const std::string& title() const { return m_title; }

 private:
  std::vector<std::string> m_node_names;
  std::map<std::string, int> m_node_index;
  std::vector<Component> m_components;
  std::vector<Component> m_vsrcs;
  std::vector<Component> m_inductors;

  std::map<std::string, DiodeModel> m_diode_models;
  std::map<std::string, MosfetModel> m_mosfet_models;
  std::map<std::string, BjtModel> m_bjt_models;

  std::map<std::string, double> m_ic;
  std::string m_ground;
  std::string m_title;
};

}  // namespace ahkab
