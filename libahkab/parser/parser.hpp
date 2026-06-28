#pragma once

#include <string>
#include <vector>

#include "core/circuit.hpp"
#include "core/types.hpp"

namespace ahkab {

struct ParsedNetlist {
  std::string title;
  Circuit circuit;
  std::vector<AnalysisStmt> analyses;
  SimOptions options;
};

ParsedNetlist parse_netlist(const std::string& content, const std::string& filename = "");

}  // namespace ahkab
