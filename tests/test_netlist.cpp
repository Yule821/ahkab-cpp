#include <gtest/gtest.h>

#include "parser/parser.hpp"

TEST(NetlistTest, ParseResistor) {
  std::string netlist =
      "Simple RC\n"
      "R1 1 0 1k\n"
      ".END\n";
  auto parsed = ahkab::parse_netlist(netlist);
  EXPECT_EQ(parsed.circuit.node_count(), 2);
  EXPECT_EQ(parsed.circuit.components().size(), 1);
}

TEST(NetlistTest, ParseWithOp) {
  std::string netlist =
      "Simple\n"
      "R1 1 0 1k\n"
      "V1 1 0 DC 5\n"
      ".OP\n"
      ".END\n";
  auto parsed = ahkab::parse_netlist(netlist);
  EXPECT_EQ(parsed.analyses.size(), 1);
  EXPECT_EQ(parsed.analyses[0].atype, ahkab::AnalysisStmt::OP);
}

TEST(NetlistTest, ParseTransient) {
  std::string netlist =
      "RC Transient\n"
      "R1 1 2 1k\n"
      "C1 2 0 1u\n"
      "V1 1 0 DC 5\n"
      ".TRAN 1u 1m\n"
      ".END\n";
  auto parsed = ahkab::parse_netlist(netlist);
  EXPECT_EQ(parsed.analyses.size(), 1);
  EXPECT_EQ(parsed.analyses[0].atype, ahkab::AnalysisStmt::TRAN);
  EXPECT_DOUBLE_EQ(parsed.analyses[0].tstep, 1e-6);
  EXPECT_DOUBLE_EQ(parsed.analyses[0].tstop, 1e-3);
}

TEST(NetlistTest, ParseDiode) {
  std::string netlist =
      "Diode Test\n"
      "D1 1 0 DMOD\n"
      ".MODEL DMOD D\n"
      ".END\n";
  auto parsed = ahkab::parse_netlist(netlist);
  EXPECT_EQ(parsed.circuit.components().size(), 1);
}
