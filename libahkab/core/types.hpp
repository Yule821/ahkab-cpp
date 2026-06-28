#pragma once

#include <complex>
#include <map>
#include <string>
#include <vector>

namespace ahkab {

enum class ComponentType {
  RESISTOR,
  CAPACITOR,
  INDUCTOR,
  VSOURCE,
  ISOURCE,
  VCVS,
  VCCS,
  CCVS,
  CCCS,
  DIODE,
  MOSFET_N,
  MOSFET_P,
  BJT_NPN,
  BJT_PNP,
};

enum class IntegrationMethod {
  IMPLICIT_EULER,
  TRAPEZOIDAL,
  GEAR2,
  GEAR3,
  GEAR4,
  GEAR5,
};

enum class SweepType {
  LINEAR,
  LOGARITHMIC,
};

struct DiodeModel {
  std::string name;
  double is = 1e-14;   // saturation current
  double n = 1.0;      // ideality factor
  double rs = 0.0;     // series resistance
  double cjo = 0.0;    // zero-bias junction capacitance
  double vj = 0.7;     // junction potential
  double m = 0.5;      // grading coefficient
  double tt = 0.0;     // transit time
  double bv = 0.0;     // breakdown voltage
  double ibv = 0.001;  // breakdown current
};

struct MosfetModel {
  std::string name;
  double vto = 0.0;       // threshold voltage
  double kp = 2e-5;       // transconductance parameter
  double gamma = 0.0;     // body effect
  double phi = 0.6;       // surface potential
  double lambda = 0.0;    // channel-length modulation
  double rd = 0.0;        // drain ohmic resistance
  double rs = 0.0;        // source ohmic resistance
  double cbd = 0.0;       // zero-bias B-D junction capacitance
  double cbs = 0.0;       // zero-bias B-S junction capacitance
  double cgso = 0.0;      // gate-source overlap capacitance
  double cgdo = 0.0;      // gate-drain overlap capacitance
  double cgbo = 0.0;      // gate-bulk overlap capacitance
  double tox = 1e-7;      // oxide thickness
  double is = 1e-14;      // bulk junction saturation current
  double pb = 0.8;        // bulk junction potential
  double cj = 0.0;        // zero-bias bulk junction capacitance
  double mj = 0.5;        // bulk junction grading coefficient
  double cjsw = 0.0;      // zero-bias sidewall capacitance
  double mjsw = 0.33;     // sidewall grading coefficient
};

struct BjtModel {
  std::string name;
  double is = 1e-16;   // transport saturation current
  double bf = 100.0;   // forward beta
  double br = 1.0;     // reverse beta
  double nf = 1.0;     // forward ideality factor
  double nr = 1.0;     // reverse ideality factor
  double vaf = 0.0;    // forward Early voltage
  double var = 0.0;    // reverse Early voltage
  double ike = 0.0;    // knee current (high-level injection)
  double ise = 0.0;    // B-E leakage saturation current
  double ne = 1.5;     // B-E leakage emission coefficient
  double isc = 0.0;    // B-C leakage saturation current
  double nc = 2.0;     // B-C leakage emission coefficient
  double rb = 0.0;     // base resistance
  double re = 0.0;     // emitter resistance
  double rc = 0.0;     // collector resistance
  double cje = 0.0;    // zero-bias B-E capacitance
  double vje = 0.75;   // B-E built-in potential
  double mje = 0.33;   // B-E grading coefficient
  double cjc = 0.0;    // zero-bias B-C capacitance
  double vjc = 0.75;   // B-C built-in potential
  double mjc = 0.33;   // B-C grading coefficient
  double tf = 0.0;     // forward transit time
  double tr = 0.0;     // reverse transit time
};

struct Component {
  std::string id;
  ComponentType type;
  std::vector<std::string> node_pins;
  std::map<std::string, double> params;
  std::string model_name;
  double dc_value = 0.0;
  double ac_magnitude = 0.0;
  double ac_phase = 0.0;
};

struct NetlistStmt {
  enum Type { COMPONENT, MODEL, ANALYSIS, OPTION, COMMENT };
  Type stmt_type;
  Component component;
  std::string raw_line;
  int line_number = 0;
};

struct AnalysisStmt {
  enum Type { OP, DC, TRAN, AC, PZ };
  Type atype;
  std::string source1;
  double start = 0.0;
  double stop = 0.0;
  double step = 0.0;
  std::string source2;
  double start2 = 0.0;
  double stop2 = 0.0;
  double step2 = 0.0;
  std::string sweep_type = "LIN";
  int npts = 0;
  double tstep = 0.0;
  double tstop = 0.0;
  double tstart = 0.0;
  double tmax = 0.0;
  std::string input_node;
  std::string output_node;
};

struct SimOptions {
  double vntol = 1e-6;
  double abstol = 1e-12;
  double reltol = 1e-3;
  int max_iter = 300;
  double gmin = 1e-12;
  int gmin_steps = 10;
  double lte_abstol = 1e-6;
  double lte_reltol = 1e-3;
  std::string integ_method = "TRAP";
  int gear_order = 2;
};

using complex_t = std::complex<double>;

constexpr double k_boltzmann = 1.380649e-23;
constexpr double q_electron = 1.602176634e-19;
constexpr double T_nominal = 300.15;
constexpr double Vt_thermal() { return k_boltzmann * T_nominal / q_electron; }

}  // namespace ahkab
