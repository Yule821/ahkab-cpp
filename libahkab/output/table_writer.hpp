#pragma once

#include <ostream>

#include "results/results.hpp"

namespace ahkab {

void write_table(const OpResult& result, std::ostream& os);

}  // namespace ahkab
