#ifndef BENIAMINO_LOSSES_H
#define BENIAMINO_LOSSES_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "numeric.h"
#include "units.h"

namespace beniamino {

template <typename T>
class LossesTable {
 public:
  LossesTable(const std::string &filename) : filename_(filename) {}

  bool loadTable() {
    std::ifstream file(filename_);
    if (!file.is_open()) {
      std::cerr << "Error: Failed to open file " << filename_ << std::endl;
      return false;
    }

    std::string line;
    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<T> row;
      std::istringstream iss(line);
      std::string cell;

      while (std::getline(iss, cell, '\t')) {
        try {
          T value = utils::convertFromString<T>(cell);
          row.push_back(value);
        } catch (const std::exception &e) {
          std::cerr << "Error: Invalid value in file " << filename_ << std::endl;
          return false;
        }
      }

      if (row.size() != 4) {
        std::cerr << "Error: Each row in the table must have 4 columns in file " << filename_
                  << std::endl;
        return false;
      }

      logE_.push_back(row[0]);
      logbeta_.push_back(row[1]);
    }

    file.close();
    return true;
  }

  T beta(T E) {
    auto logE = std::log10(E / SI::eV);
    double logbeta = 0;
    if (logE < logE_.front())
      logbeta = logbeta_.front();
    else if (logE > logE_.back())
      logbeta = logbeta_.back();
    else
      logbeta = utils::interpolate(logE, logE_, logbeta_);
    return std::pow(10., logbeta) / SI::year;
  }

  T dbdE(T E) {
    auto logE = std::log10(E / SI::eV);
    if (logE <= logE_.front() || logE >= logE_.back()) return beta(E);
    size_t i = std::lower_bound(logE_.begin(), logE_.end(), logE) - logE_.begin();
    assert(logE <= logE_.at(i) && logE >= logE_.at(i - 1));
    return beta(E) * (1. + (logbeta_[i] - logbeta_[i - 1]) / (logE_[i] - logE_[i - 1]));
  }

 private:
  std::string filename_;
  std::vector<T> logE_;
  std::vector<T> logbeta_;
};

}  // namespace beniamino

#endif