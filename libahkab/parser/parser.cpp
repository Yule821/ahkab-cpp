#include "parser/parser.hpp"

#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>

#include "core/errors.hpp"
#include "parser/lexer.hpp"

namespace ahkab {

static double parse_value(const std::string& s) {
  std::string clean;
  for (char c : s) {
    if ((c >= '0' && c <= '9') || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E') {
      clean += c;
    }
  }
  if (s.find("MEG") != std::string::npos) return std::stod(clean) * 1e6;
  if (s.find("MIL") != std::string::npos) return std::stod(clean) * 25.4e-6;
  double v = std::stod(clean);
  char last = '\0';
  for (char c : s) {
    if ((c < '0' || c > '9') && c != '.' && c != '-' && c != '+' && c != 'e' && c != 'E') {
      last = c;
    }
  }
  if (last == '\0') return v;
  switch (std::toupper(last)) {
    case 'T': return v * 1e12;
    case 'G': return v * 1e9;
    case 'K': return v * 1e3;
    case 'M': return v * 1e-3;
    case 'U': return v * 1e-6;
    case 'N': return v * 1e-9;
    case 'P': return v * 1e-12;
    case 'F': return v * 1e-15;
    case 'H': return v;
    case 'R': return v;
    case 'L': return v;
    case 'V': return v;
    case 'A': return v;
    case 'S': return v;
    default: return v;
  }
}

ParsedNetlist parse_netlist(const std::string& content, const std::string& filename) {
  ParsedNetlist result;
  auto tokens = tokenize(content, filename);

  if (tokens.empty()) {
    throw NetlistError("Empty netlist", filename, 1);
  }

  size_t idx;
  if (!tokens[0].value.empty() && tokens[0].value[0] != '.') {
    result.title = "";
    int title_line = tokens[0].line;
    size_t t = 0;
    while (t < tokens.size() && tokens[t].line == title_line) {
      if (!result.title.empty()) result.title += " ";
      result.title += tokens[t].value;
      t++;
    }
    idx = t;
  } else {
    result.title = tokens[0].value;
    idx = 1;
  }
  while (idx < tokens.size()) {
    std::string first = tokens[idx].value;

    if (first == ".MODEL") {
      if (idx + 3 >= tokens.size()) {
        throw NetlistError("Incomplete .MODEL statement", filename,
                           tokens[idx].line);
      }
      std::string model_name = tokens[++idx].value;
      std::string model_type = tokens[++idx].value;
      std::string params_str;
      idx++;
      while (idx < tokens.size() && tokens[idx].value[0] != '.') {
        if (!params_str.empty()) params_str += " ";
        params_str += tokens[idx].value;
        idx++;
      }

      if (model_type == "D") {
        DiodeModel dm;
        dm.name = model_name;
        result.circuit.add_diode_model(model_name, dm);
      } else if (model_type == "NMOS" || model_type == "PMOS") {
        MosfetModel mm;
        mm.name = model_name;
        result.circuit.add_mosfet_model(model_name, mm);
      } else if (model_type == "NPN" || model_type == "PNP") {
        BjtModel bm;
        bm.name = model_name;
        result.circuit.add_bjt_model(model_name, bm);
      }
      continue;
    }

    if (first == ".OP") {
      AnalysisStmt a;
      a.atype = AnalysisStmt::OP;
      result.analyses.push_back(a);
      idx++;
      continue;
    }

    if (first == ".DC") {
      if (idx + 4 >= tokens.size()) {
        throw NetlistError("Incomplete .DC statement", filename,
                           tokens[idx].line);
      }
      AnalysisStmt a;
      a.atype = AnalysisStmt::DC;
      a.source1 = tokens[++idx].value;
      a.start = parse_value(tokens[++idx].value);
      a.stop = parse_value(tokens[++idx].value);
      a.step = parse_value(tokens[++idx].value);
      idx++;
      result.analyses.push_back(a);
      continue;
    }

    if (first == ".TRAN") {
      if (idx + 3 >= tokens.size()) {
        throw NetlistError("Incomplete .TRAN statement", filename,
                           tokens[idx].line);
      }
      AnalysisStmt a;
      a.atype = AnalysisStmt::TRAN;
      a.tstep = parse_value(tokens[++idx].value);
      a.tstop = parse_value(tokens[++idx].value);
      a.tstart = 0.0;
      idx++;
      if (idx < tokens.size() && tokens[idx].value[0] != '.') {
        a.tstart = std::stod(tokens[idx].value);
        idx++;
      }
      result.analyses.push_back(a);
      continue;
    }

    if (first == ".AC") {
      if (idx + 4 >= tokens.size()) {
        throw NetlistError("Incomplete .AC statement", filename,
                           tokens[idx].line);
      }
      AnalysisStmt a;
      a.atype = AnalysisStmt::AC;
      a.sweep_type = tokens[++idx].value;
      a.npts = std::stoi(tokens[++idx].value);
      a.start = parse_value(tokens[++idx].value);
      a.stop = parse_value(tokens[++idx].value);
      idx++;
      result.analyses.push_back(a);
      continue;
    }

    if (first == ".PZ") {
      if (idx + 3 >= tokens.size()) {
        throw NetlistError("Incomplete .PZ statement", filename,
                           tokens[idx].line);
      }
      AnalysisStmt a;
      a.atype = AnalysisStmt::PZ;
      a.input_node = tokens[++idx].value;
      a.output_node = tokens[++idx].value;
      idx++;
      result.analyses.push_back(a);
      continue;
    }

    if (first == ".END") {
      break;
    }

    if (first == ".TITLE") {
      if (idx + 1 < tokens.size()) {
        result.title = tokens[++idx].value;
      }
      idx++;
      continue;
    }

    char ctype = std::toupper(first[0]);
    std::string id = first;
    idx++;

    switch (ctype) {
      case 'R': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        double val = parse_value(tokens[idx + 2].value);
        result.circuit.add_resistor(id, n1, n2, val);
        idx += 3;
        break;
      }
      case 'C': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        double val = parse_value(tokens[idx + 2].value);
        result.circuit.add_capacitor(id, n1, n2, val);
        idx += 3;
        break;
      }
      case 'L': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        double val = parse_value(tokens[idx + 2].value);
        result.circuit.add_inductor(id, n1, n2, val);
        idx += 3;
        break;
      }
      case 'V': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        double dc = 0.0;
        double ac_mag = 0.0, ac_phase = 0.0;
        idx += 2;
        if (idx < tokens.size() && tokens[idx].value == "DC") {
          idx++;
          if (idx < tokens.size()) {
            dc = std::stod(tokens[idx].value);
            idx++;
          }
          if (idx < tokens.size() && tokens[idx].value == "AC") {
            idx++;
            if (idx < tokens.size()) ac_mag = std::stod(tokens[idx].value);
            idx++;
            if (idx < tokens.size()) ac_phase = std::stod(tokens[idx].value);
            idx++;
          }
        } else if (idx < tokens.size() && tokens[idx].value == "AC") {
          idx++;
          if (idx < tokens.size()) ac_mag = std::stod(tokens[idx].value);
          idx++;
          if (idx < tokens.size()) ac_phase = std::stod(tokens[idx].value);
          idx++;
        } else {
          dc = std::stod(tokens[idx].value);
          idx++;
          if (idx < tokens.size() && tokens[idx].value == "AC") {
            idx++;
            if (idx < tokens.size()) ac_mag = std::stod(tokens[idx].value);
            idx++;
            if (idx < tokens.size()) ac_phase = std::stod(tokens[idx].value);
            idx++;
          }
        }
        result.circuit.add_vsource(id, n1, n2, dc, ac_mag, ac_phase);
        break;
      }
      case 'I': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        double dc = 0.0;
        idx += 2;
        if (idx < tokens.size() && tokens[idx].value == "DC") {
          idx++;
          if (idx < tokens.size()) dc = std::stod(tokens[idx].value);
          idx++;
        } else {
          dc = std::stod(tokens[idx].value);
          idx++;
        }
        result.circuit.add_isource(id, n1, n2, dc);
        break;
      }
      case 'D': {
        if (idx + 2 >= tokens.size()) break;
        std::string anode = tokens[idx].value;
        std::string cathode = tokens[idx + 1].value;
        std::string model = tokens[idx + 2].value;
        result.circuit.add_diode(id, anode, cathode, model);
        idx += 3;
        break;
      }
      case 'M': {
        if (idx + 4 >= tokens.size()) break;
        std::string drain = tokens[idx].value;
        std::string gate = tokens[idx + 1].value;
        std::string source = tokens[idx + 2].value;
        std::string bulk = tokens[idx + 3].value;
        std::string model = tokens[idx + 4].value;
        result.circuit.add_mosfet(id, drain, gate, source, bulk, model,
                                  ComponentType::MOSFET_N);
        idx += 5;
        break;
      }
      case 'Q': {
        if (idx + 3 >= tokens.size()) break;
        std::string collector = tokens[idx].value;
        std::string base = tokens[idx + 1].value;
        std::string emitter = tokens[idx + 2].value;
        std::string model = tokens[idx + 3].value;
        result.circuit.add_bjt(id, collector, base, emitter, model,
                               ComponentType::BJT_NPN);
        idx += 4;
        break;
      }
      case 'E': {
        if (idx + 3 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        std::string cn1 = tokens[idx + 2].value;
        std::string cn2 = tokens[idx + 3].value;
        double gain = std::stod(tokens[idx + 4].value);
        result.circuit.add_vcvs(id, n1, n2, cn1, cn2, gain);
        idx += 5;
        break;
      }
      case 'F': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        std::string vsource = tokens[idx + 2].value;
        double gain = std::stod(tokens[idx + 3].value);
        result.circuit.add_cccs(id, n1, n2, vsource, gain);
        idx += 4;
        break;
      }
      case 'G': {
        if (idx + 3 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        std::string cn1 = tokens[idx + 2].value;
        std::string cn2 = tokens[idx + 3].value;
        double gain = std::stod(tokens[idx + 4].value);
        result.circuit.add_vccs(id, n1, n2, cn1, cn2, gain);
        idx += 5;
        break;
      }
      case 'H': {
        if (idx + 2 >= tokens.size()) break;
        std::string n1 = tokens[idx].value;
        std::string n2 = tokens[idx + 1].value;
        std::string vsource = tokens[idx + 2].value;
        double gain = std::stod(tokens[idx + 3].value);
        result.circuit.add_ccvs(id, n1, n2, vsource, gain);
        idx += 4;
        break;
      }
      default:
        idx++;
        break;
    }
  }

  return result;
}

}  // namespace ahkab
