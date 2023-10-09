#include "losses.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "simprop/core/units.h"
#include "simprop/utils/numeric.h"
#include "utils.h"

namespace beniamino {

bool LossesTable::loadTable(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Failed to open file " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<double> row;
    std::istringstream iss(line);
    std::string cell;

    while (std::getline(iss, cell, '\t')) {
      try {
        double value = utils::convertFromString<double>(cell);
        row.push_back(value);
      } catch (const std::exception& e) {
        std::cerr << "Error: Invalid value in file " << filename << std::endl;
        return false;
      }
    }

    if (row.size() != 3) {
      std::cerr << "Error: Each row in the table must have 4 columns in file " << filename
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

double LossesTable::beta(double E) const {
  auto logE = std::log10(E / SI::eV);
  double logbeta = 0;
  if (logE < logE_.front())
    logbeta = logbeta_.front();
  else if (logE > logE_.back())
    logbeta = logbeta_.back();
  else
    logbeta = simprop::utils::interpolate(logE, logE_, logbeta_);
  return std::pow(10., logbeta) / SI::year;
}

double LossesTable::dbdE(double E) const {
  auto logE = std::log10(E / SI::eV);
  double logdbde = 0;
  if (logE < logE_.front())
    logdbde = logdbdE_.front();
  else if (logE > logE_.back())
    logdbde = logdbdE_.back();
  else
    logdbde = simprop::utils::interpolate(logE, logE_, logdbdE_);
  return std::pow(10., logdbde) / SI::year;
}

}  // namespace beniamino