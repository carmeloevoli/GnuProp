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

      if (row.size() != 3) {
        std::cerr << "Error: Each row in the table must have 4 columns in file " << filename_
                  << std::endl;
        return false;
      }

      logE_.push_back(row[0]);
      logbeta_.push_back(row[1]);
      logdbdE_.push_back(row[2]);
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
    double logdbde = 0;
    if (logE < logE_.front())
      logdbde = logdbdE_.front();
    else if (logE > logE_.back())
      logdbde = logdbdE_.back();
    else
      logdbde = utils::interpolate(logE, logE_, logdbdE_);
    return std::pow(10., logdbde) / SI::year;
  }

 private:
  std::string filename_;
  std::vector<T> logE_;
  std::vector<T> logbeta_;
  std::vector<T> logdbdE_;
};

}  // namespace beniamino

#endif