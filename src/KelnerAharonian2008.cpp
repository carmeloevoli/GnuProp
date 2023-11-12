#include "KelnerAharonian2008.h"

#include <fstream>
#include <iostream>

#include "simprop/utils/numeric.h"
#include "utils.h"

namespace KelnerAharonian2008 {

bool SecondarySpectrum::loadTable(const std::string& filename) {
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

    while (std::getline(iss, cell, ',')) {
      try {
        double value = utils::convertFromString<double>(cell);
        row.push_back(value);
      } catch (const std::exception& e) {
        std::cerr << "Error: Invalid value in file " << filename << std::endl;
        return false;
      }
    }

    if (row.size() != 4) {
      std::cerr << "Error: Each row in the table must have 4 columns in file " << filename
                << std::endl;
      return false;
    }

    m_rho_table.push_back(row[0]);
    m_s_table.push_back(row[1]);
    m_delta_table.push_back(row[2]);
    m_lnB_table.push_back(std::log(row[3] * SI::cm3 / SI::sec));
  }

  file.close();
  return true;
}

double SecondarySpectrum::B(double rho) const {
  double value = 0;
  if (rho > m_rho_table.front() && rho < m_rho_table.back()) {
    auto lnB = simprop::utils::cspline(rho, m_rho_table, m_lnB_table);
    value = std::exp(lnB);
  }
  return value;
}

double SecondarySpectrum::s(double rho) const {
  double value = 0;
  if (rho > m_rho_table.front() && rho < m_rho_table.back()) {
    return simprop::utils::cspline(rho, m_rho_table, m_s_table);
  }
  return value;
}

double SecondarySpectrum::delta(double rho) const {
  double value = 0;
  if (rho > m_rho_table.front() && rho < m_rho_table.back()) {
    return simprop::utils::cspline(rho, m_rho_table, m_delta_table);
  }
  return value;
}

double SecondarySpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;
  const auto rho = eta / m_eta_0;
  const auto _B = B(rho);
  const auto _s = s(rho);
  const auto _delta = delta(rho);
  const auto _psi = psi(rho);

  const auto _xPrimeMinus = xPrimeMinus(eta);
  const auto _xPrimePlus = xPrimePlus(eta);

  const auto yPrime = (x - _xPrimeMinus) / (_xPrimePlus - _xPrimeMinus);

  double value = 0.;
  if (x < _xPrimeMinus) {
    value = _B * std::pow(M_LN2, _psi);
  } else if (x < _xPrimePlus) {
    value = _B;
    value *= std::exp(-_s * std::pow(std::log(x / _xPrimeMinus), _delta));
    value *= std::pow(std::log(2. / (1. + pow2(yPrime))), _psi);
  }
  return value;
}

double SecondarySpectrum::xMinus(double eta) const {
  auto value = 1. / 2. / (1. + eta);
  value *= eta + m_r2 - std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
  return value;
}

double SecondarySpectrum::xPlus(double eta) const {
  auto value = 1. / 2. / (1. + eta);
  value *= eta + m_r2 + std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
  return value;
}

AntiNuMuSpectrum::AntiNuMuSpectrum() { loadTable("data/xsecs_KA2018_antiNuMu.txt"); }

double AntiNuMuSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double AntiNuMuSpectrum::xPrimeMinus(double eta) const { return xMinus(eta) / 4.; }

double AntiNuMuSpectrum::xPrimePlus(double eta) const { return xPlus(eta); }

NuMuSpectrum::NuMuSpectrum() { loadTable("data/xsecs_KA2018_nuMu.txt"); }

double NuMuSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double NuMuSpectrum::xPrimeMinus(double eta) const { return 0.427 * xMinus(eta); }

double NuMuSpectrum::xPrimePlus(double eta) const {
  const double rho = eta / m_eta_0;
  if (rho < 2.14) {
    return 0.427 * xPlus(eta);
  } else if (rho < 10.) {
    return (0.427 + 0.0729 * (rho - 2.14)) * xPlus(eta);
  } else {
    return xPlus(eta);
  }
}

NuElectronSpectrum::NuElectronSpectrum() { loadTable("data/xsecs_KA2018_nuElectron.txt"); }

double NuElectronSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double NuElectronSpectrum::xPrimeMinus(double eta) const { return xMinus(eta) / 4.; }

double NuElectronSpectrum::xPrimePlus(double eta) const { return xPlus(eta); }

AntiNuElectronSpectrum::AntiNuElectronSpectrum() {
  loadTable("data/xsecs_KA2018_antiNuElectron.txt");
}

double AntiNuElectronSpectrum::psi(double rho) const {
  return (rho > 4.) ? 6. * (1. - std::exp(1.5 * (4. - rho))) : 0.;
}

double AntiNuElectronSpectrum::xPrimeMinus(double eta) const {
  double value = 1. / 2. / (1. + eta);
  value *= eta - 2. * m_r - std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)));
  return value;
}

double AntiNuElectronSpectrum::xPrimePlus(double eta) const {
  double value = 1. / 2. / (1. + eta);
  value *= eta - 2. * m_r + std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)));
  return value;
}

}  // namespace KelnerAharonian2008
