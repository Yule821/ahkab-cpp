#include <cmath>
#include <gtest/gtest.h>

#include "analysis/op_analysis.hpp"
#include "core/circuit.hpp"

TEST(OpAnalysisTest, ResistorDivider) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_vsource("V1", "1", "0", 5.0);
  c.add_resistor("R1", "1", "2", 1000.0);
  c.add_resistor("R2", "2", "0", 1000.0);

  ahkab::NrConfig cfg;
  cfg.max_iter = 100;
  auto result = ahkab::run_op(c, cfg);

  EXPECT_TRUE(result.meta.converged);

  EXPECT_NEAR(result.node_voltages["V(1)"], 5.0, 1e-3);
  EXPECT_NEAR(result.node_voltages["V(2)"], 2.5, 1e-3);
}

TEST(OpAnalysisTest, VoltageDividerWithDiode) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_vsource("V1", "1", "0", 10.0);
  c.add_resistor("R1", "1", "2", 1000.0);

  ahkab::DiodeModel dm;
  dm.name = "DMOD";
  dm.is = 1e-14;
  dm.n = 1.0;
  c.add_diode_model("DMOD", dm);
  c.add_diode("D1", "2", "0", "DMOD");

  ahkab::NrConfig cfg;
  cfg.max_iter = 200;
  auto result = ahkab::run_op(c, cfg);

  EXPECT_TRUE(result.meta.converged);
  EXPECT_NEAR(result.node_voltages["V(1)"], 10.0, 1e-2);
  EXPECT_GT(result.node_voltages["V(2)"], 0.5);
}
