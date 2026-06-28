#pragma once

#include <stdexcept>
#include <string>

namespace ahkab {

class Error : public std::runtime_error {
 public:
  explicit Error(const std::string& msg) : std::runtime_error(msg) {}
};

class NetlistError : public Error {
 public:
  NetlistError(const std::string& msg, const std::string& file,
               int line)
      : Error(file + ":" + std::to_string(line) + ": " + msg),
        filename(file),
        line_number(line) {}

  std::string filename;
  int line_number;
};

class ConvergenceError : public Error {
 public:
  ConvergenceError(const std::string& msg, double residual, int iterations)
      : Error(msg), residual_norm(residual), iter_count(iterations) {}

  double residual_norm;
  int iter_count;
};

class TimeStepError : public Error {
 public:
  TimeStepError(const std::string& msg, double t)
      : Error("at t=" + std::to_string(t) + ": " + msg), current_time(t) {}

  double current_time;
};

class SingularMatrixError : public Error {
 public:
  explicit SingularMatrixError(const std::string& msg) : Error(msg) {}
};

}  // namespace ahkab
