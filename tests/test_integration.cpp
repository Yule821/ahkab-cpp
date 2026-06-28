#include <cmath>
#include <gtest/gtest.h>

#include "analysis/ac_analysis.hpp"
#include "analysis/dc_analysis.hpp"
#include "analysis/op_analysis.hpp"

TEST(IntegrationTest, ResistorDividerFullFlow) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_vsource("V1", "1", "0", 10.0);
  c.add_resistor("R1", "1", "2", 1000.0);
  c.add_resistor("R2", "2", "0", 2000.0);

  ahkab::NrConfig cfg;
  cfg.max_iter = 100;
  auto op = ahkab::run_op(c, cfg);
  EXPECT_TRUE(op.meta.converged);
  EXPECT_NEAR(op.node_voltages["V(2)"], 10.0 * 2000.0 / 3000.0, 1e-2);
}

TEST(IntegrationTest, DcSweep) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_vsource("V1", "1", "0", 1.0);
  c.add_resistor("R1", "1", "2", 1000.0);
  c.add_resistor("R2", "2", "0", 1000.0);

  auto result = ahkab::run_dc(c, "V1", 1.0, 6.0, 1.0);
  EXPECT_GE(result.sweep_values.size(), 5);
  EXPECT_NEAR(result.sweep_values[0], 1.0, 1e-9);
  EXPECT_NEAR(result.sweep_values.back(), 6.0, 1e-9);
}

TEST(IntegrationTest, AcAnalysis) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_vsource("V1", "1", "0", 1.0, 1.0);
  c.add_resistor("R1", "1", "2", 1000.0);
  c.add_capacitor("C1", "2", "0", 1e-6);

  ahkab::AcConfig acfg;
  acfg.fstart = 10.0;
  acfg.fstop = 1000.0;
  acfg.npts = 5;
  acfg.sweep_type = ahkab::SweepType::LOGARITHMIC;

  auto result = ahkab::run_ac(c, acfg);
  EXPECT_EQ(result.frequencies.size(), 5);
  EXPECT_GT(std::abs(result.values[0][1]), 0.0);
}
