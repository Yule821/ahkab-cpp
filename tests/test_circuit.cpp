#include <gtest/gtest.h>

#include "core/circuit.hpp"

TEST(CircuitTest, AddResistor) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_resistor("R1", "1", "0", 1000.0);
  EXPECT_EQ(c.node_count(), 2);
  EXPECT_EQ(c.mna_dimension(), 1);
}

TEST(CircuitTest, AddVsource) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_vsource("V1", "1", "0", 5.0);
  EXPECT_EQ(c.num_vsrcs(), 1);
  EXPECT_EQ(c.mna_dimension(), 2);
}

TEST(CircuitTest, ResistorDivider) {
  ahkab::Circuit c;
  c.set_ground("0");
  c.add_resistor("R1", "1", "2", 1000.0);
  c.add_resistor("R2", "2", "0", 1000.0);
  EXPECT_EQ(c.node_count(), 3);
}
