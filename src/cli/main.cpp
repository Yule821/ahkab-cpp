#include <fstream>
#include <iostream>
#include <string>

#include "analysis/ac_analysis.hpp"
#include "analysis/dc_analysis.hpp"
#include "analysis/op_analysis.hpp"
#include "analysis/pz_analysis.hpp"
#include "analysis/transient_analysis.hpp"
#include "core/errors.hpp"
#include "output/csv_writer.hpp"
#include "parser/parser.hpp"

static void print_usage() {
  std::cerr << "Usage: ahkab_cli [-o output] <netlist_file>\n";
}

int main(int argc, char* argv[]) {
  std::string output_file;
  std::string input_file;

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-o" && i + 1 < argc) {
      output_file = argv[++i];
    } else if (arg == "-h" || arg == "--help") {
      print_usage();
      return 0;
    } else {
      input_file = arg;
    }
  }

  if (input_file.empty()) {
    print_usage();
    return 1;
  }

  try {
    std::ifstream ifs(input_file);
    if (!ifs) {
      std::cerr << "Error: cannot open " << input_file << "\n";
      return 1;
    }
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    auto parsed = ahkab::parse_netlist(content, input_file);
    auto& circuit = parsed.circuit;
    circuit.set_ground("0");

    std::ostream* out = &std::cout;
    std::ofstream ofs;
    if (!output_file.empty()) {
      ofs.open(output_file);
      out = &ofs;
    }

    for (const auto& analysis : parsed.analyses) {
      switch (analysis.atype) {
        case ahkab::AnalysisStmt::OP: {
          auto result = ahkab::run_op(circuit);
          ahkab::write_table(result, *out);
          break;
        }
        case ahkab::AnalysisStmt::DC: {
          auto result = ahkab::run_dc(circuit, analysis.source1,
                                       analysis.start, analysis.stop,
                                       analysis.step);
          *out << "DC Sweep: " << analysis.source1
               << " from " << analysis.start
               << " to " << analysis.stop
               << " step " << analysis.step << "\n";
          for (size_t i = 0; i < result.sweep_values.size(); ++i) {
            *out << result.sweep_values[i];
            for (double v : result.data[i]) {
              *out << " " << v;
            }
            *out << "\n";
          }
          break;
        }
        case ahkab::AnalysisStmt::TRAN: {
          ahkab::TranConfig tcfg;
          tcfg.tstart = analysis.tstart;
          tcfg.tstop = analysis.tstop;
          tcfg.tstep = analysis.tstep;
          auto result = ahkab::run_tran(circuit, tcfg);
          ahkab::write_csv(result, *out);
          break;
        }
        case ahkab::AnalysisStmt::AC: {
          ahkab::AcConfig acfg;
          acfg.fstart = analysis.start;
          acfg.fstop = analysis.stop;
          acfg.npts = analysis.npts;
          acfg.sweep_type = (analysis.sweep_type == "LIN")
                                 ? ahkab::SweepType::LINEAR
                                 : ahkab::SweepType::LOGARITHMIC;
          auto result = ahkab::run_ac(circuit, acfg);
          *out << "freq";
          for (const auto& name : result.variable_names) {
            *out << ",mag(" << name << "),phase(" << name << ")";
          }
          *out << "\n";
          for (size_t i = 0; i < result.frequencies.size(); ++i) {
            *out << result.frequencies[i];
            for (const auto& v : result.values[i]) {
              *out << "," << std::abs(v) << "," << std::arg(v);
            }
            *out << "\n";
          }
          break;
        }
        case ahkab::AnalysisStmt::PZ: {
          ahkab::PzConfig pcfg;
          pcfg.input_source = analysis.input_node;
          pcfg.output_node = analysis.output_node;
          auto result = ahkab::run_pz(circuit, pcfg);
          *out << "Poles:\n";
          for (const auto& p : result.poles) {
            *out << p.real() << " + j" << p.imag() << "\n";
          }
          *out << "Zeros:\n";
          for (const auto& z : result.zeros) {
            *out << z.real() << " + j" << z.imag() << "\n";
          }
          break;
        }
      }
    }

  } catch (const ahkab::Error& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
