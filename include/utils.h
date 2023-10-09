#ifndef BENIAMINO_UTILS_H
#define BENIAMINO_UTILS_H

#include <sstream>

namespace utils {

template <typename T>
T convertFromString(const std::string &str) {
  std::istringstream iss(str);
  T value;
  iss >> value;
  if (!iss.eof() || iss.fail()) {
    throw std::invalid_argument("Failed to convert string to the specified type");
  }
  return value;
}

}  // namespace utils

#endif