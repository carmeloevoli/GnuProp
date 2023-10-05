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
    table_.resize(3);

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

      table_[0].push_back(row[0]);
      table_[1].push_back(row[1]);
      table_[2].push_back(row[2]);
    }

    file.close();
    return true;
  }

  T beta(T E) {
    auto logE = std::log10(E / SI::eV);
    auto logb = utils::interpolate(logE, table_[0], table_[1]);
    return std::pow(10., logb) / SI::year;
  }

  T dbdE(T E) {
    auto logE = std::log10(E / SI::eV);
    auto logb = utils::interpolate(logE, table_[0], table_[2]);
    return std::pow(10., logb) / SI::year;
  }

 private:
  std::string filename_;
  std::vector<std::vector<T>> table_;
};

}  // namespace beniamino

#endif