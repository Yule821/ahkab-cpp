#pragma once

#include <string>
#include <vector>

namespace ahkab {

struct Token {
  std::string value;
  int line;
};

std::vector<Token> tokenize(const std::string& content, const std::string& filename);

}  // namespace ahkab
