# ahkab-cpp

C++ SPICE-style circuit simulator, ported from [ahkab](https://github.com/ahkab/ahkab) (Python).

## Build

```bash
mkdir build && cd build
cmake .. -DAHKAB_BUILD_TESTS=ON
make -j$(nproc)
```

Dependencies: Eigen3, GTest (optional, for tests).

## Usage

### CLI

```bash
./build/src/cli/ahkab_cli test_rc.cir
```

### C++ API

```cpp
#include "core/circuit.hpp"
#include "analysis/op_analysis.hpp"

ahkab::Circuit c;
c.set_ground("0");
c.add_vsource("V1", "1", "0", 10.0);
c.add_resistor("R1", "1", "2", 1000.0);
c.add_resistor("R2", "2", "0", 2000.0);

auto result = ahkab::run_op(c);
std::cout << "V(2) = " << result.node_voltages.at("V(2)") << std::endl;
```

### Analysis types

| Command | Description |
|---------|-------------|
| `.OP` | DC operating point |
| `.DC Vsrc start stop step` | DC sweep |
| `.TRAN tstep tstop` | Transient analysis |
| `.AC DEC npts fstart fstop` | AC frequency sweep |
| `.PZ in_node out_node` | Pole-zero analysis |

### Supported components

R, C, L, V, I, E(VCVS), G(VCCS), H(CCVS), F(CCCS), D(diode), M(MOSFET), Q(BJT)

## Test

```bash
cd build && ctest
```
