#include "parser/lexer.hpp"

#include <cctype>

#include "core/errors.hpp"

namespace ahkab {

std::vector<Token> tokenize(const std::string& content, const std::string& filename) {
  std::vector<Token> tokens;
  int line = 1;
  size_t i = 0;

  while (i < content.size()) {
    while (i < content.size() && (content[i] == ' ' || content[i] == '\t')) {
      i++;
    }

    if (i >= content.size()) break;

    if (content[i] == '\n') {
      line++;
      i++;
      continue;
    }

    if (content[i] == '\r') {
      i++;
      continue;
    }

    if (content[i] == '*') {
      while (i < content.size() && content[i] != '\n') i++;
      continue;
    }

    if (content[i] == '+') {
      i++;
      while (i < content.size() && content[i] == ' ') i++;
      if (i < content.size() && content[i] == '\n') {
        line++;
        i++;
      }
      continue;
    }

    std::string value;
    if (content[i] == '"' || content[i] == '\'') {
      char quote = content[i++];
      while (i < content.size() && content[i] != quote) {
        value += content[i++];
      }
      if (i < content.size()) i++;
    } else {
      while (i < content.size() && !std::isspace(content[i])) {
        value += content[i++];
      }
    }

    if (!value.empty()) {
      tokens.push_back({value, line});
    }
  }

  return tokens;
}

}  // namespace ahkab
