#pragma once

#include <ostream>

#include "results/results.hpp"

namespace ahkab {

void write_csv(const TranResult& result, std::ostream& os);
void write_table(const OpResult& result, std::ostream& os);

}  // namespace ahkab
