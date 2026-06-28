#pragma once

#include <complex>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "core/types.hpp"
#include "mna/nr_solver.hpp"

namespace ahkab {

struct OpResult {
  std::map<std::string, double> node_voltages;
  std::map<std::string, double> branch_currents;
  NrResult meta;
};

struct DcResult {
  std::vector<std::string> variable_names;
  std::vector<double> sweep_values;
  std::vector<std::vector<double>> data;
};

struct TranResult {
  std::vector<double> time;
  std::vector<std::string> variable_names;
  std::vector<std::vector<double>> values;
};

struct AcResult {
  std::vector<double> frequencies;
  std::vector<std::string> variable_names;
  std::vector<std::vector<complex_t>> values;
};

struct PzResult {
  std::vector<complex_t> poles;
  std::vector<complex_t> zeros;
};

}  // namespace ahkab
