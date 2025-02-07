#include "interactions/PhotoPion.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "simprop/utils/numeric.h"
#include "utils.h"

namespace interactions {

std::vector<std::vector<double>> loadTables(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Error: Failed to open file " + filename);
  }

  std::vector<double> lnrho, s, delta, lnB;
  std::string line;
  size_t lineNumber = 0;

  while (std::getline(file, line)) {
    ++lineNumber;

    // Skip empty lines and comments
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    std::vector<double> row;
    std::string cell;

    while (std::getline(iss, cell, ',')) {
      try {
        row.push_back(utils::convertFromString<double>(cell));
      } catch (const std::exception& e) {
        throw std::runtime_error("Error in file " + filename + " at line " +
                                 std::to_string(lineNumber) + ": " + e.what());
      }
    }

    // Ensure each row has exactly 4 values
    if (row.size() != 4) {
      throw std::runtime_error("Error in file " + filename + " at line " +
                               std::to_string(lineNumber) + ": Expected 4 columns, got " +
                               std::to_string(row.size()));
    }

    // Store converted values
    lnrho.push_back(std::log(row[0]));
    s.push_back(row[1]);
    delta.push_back(row[2]);
    lnB.push_back(std::log(row[3] * SI::cm3 / SI::sec));
  }

  return {lnrho, s, delta, lnB};
}

KelnerAharonian2008::KelnerAharonian2008(const std::string& filename) {
  auto v = loadTables(filename);

  m_lnrho_table = std::move(v[0]);
  m_s_table = std::move(v[1]);
  m_delta_table = std::move(v[2]);
  m_lnB_table = std::move(v[3]);
}

double KelnerAharonian2008::B(double rho) const {
  double value = 0;
  auto ln_rho = std::log(rho);
  if (ln_rho > m_lnrho_table.front() && ln_rho < m_lnrho_table.back()) {
    auto lnB = simprop::utils::interpolate(ln_rho, m_lnrho_table, m_lnB_table);
    value = std::exp(lnB);
  }
  return value;
}

double KelnerAharonian2008::s(double rho) const {
  double value = 0;
  auto ln_rho = std::log(rho);
  if (ln_rho > m_lnrho_table.front() && ln_rho < m_lnrho_table.back()) {
    return simprop::utils::interpolate(ln_rho, m_lnrho_table, m_s_table);
  }
  return value;
}

double KelnerAharonian2008::delta(double rho) const {
  double value = 0;
  auto ln_rho = std::log(rho);
  if (ln_rho > m_lnrho_table.front() && ln_rho < m_lnrho_table.back()) {
    return simprop::utils::interpolate(ln_rho, m_lnrho_table, m_delta_table);
  }
  return value;
}

double KelnerAharonian2008::xMinus(double eta) const {  // Eq. 19
  auto value = 1. / 2. / (1. + eta);
  value *= eta + m_r2 - std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
  return value;
}

double KelnerAharonian2008::xPlus(double eta) const {  // Eq. 19
  auto value = 1. / 2. / (1. + eta);
  value *= eta + m_r2 + std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
  return value;
}

double func(double x, double xMinus, double xPlus, double s, double delta, double B, double psi) {
  auto y = (x - xMinus) / (xPlus - xMinus);  // Eq. 28
  if (x < xMinus) {                          // Eq. 29
    return B * std::pow(M_LN2, psi);
  } else if (x < xPlus) {  // Eq. 27
    return B * std::exp(-s * std::pow(std::log(x / xMinus), delta)) *
           std::pow(std::log(2. / (1. + pow2(y))), psi);
  }
  return 0.;
}

AntiNuMuSpectrum::AntiNuMuSpectrum() : KelnerAharonian2008("data/xsecs_KA2018_antiNuMu.txt") {}

double AntiNuMuSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double AntiNuMuSpectrum::xPrimeMinus(double eta) const { return xMinus(eta) / 4.; }

double AntiNuMuSpectrum::xPrimePlus(double eta) const { return xPlus(eta); }

double AntiNuMuSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xPrimeMinus = xPrimeMinus(eta);
  const double _xPrimePlus = xPrimePlus(eta);

  return func(x, _xPrimeMinus, _xPrimePlus, _s, _delta, _B, _psi);
}

NuMuSpectrum::NuMuSpectrum() : KelnerAharonian2008("data/xsecs_KA2018_nuMu.txt") {}

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

double NuMuSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xPrimeMinus = xPrimeMinus(eta);
  const double _xPrimePlus = xPrimePlus(eta);

  return func(x, _xPrimeMinus, _xPrimePlus, _s, _delta, _B, _psi);
}

AntiNuElectronSpectrum::AntiNuElectronSpectrum()
    : KelnerAharonian2008("data/xsecs_KA2018_antiNuElectron.txt") {}

double AntiNuElectronSpectrum::psi(double rho) const {
  return (rho > 4.) ? 6. * (1. - std::exp(1.5 * (4. - rho))) : 0.;
}

double AntiNuElectronSpectrum::xPrimeMinus(double eta) const {
  double value = 1. / 2. / (1. + eta);
  value *= eta - 2. * m_r - std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)));
  return value / 2.;
}

double AntiNuElectronSpectrum::xPrimePlus(double eta) const {
  double value = 1. / 2. / (1. + eta);
  value *= eta - 2. * m_r + std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)));
  return value;
}

double AntiNuElectronSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xPrimeMinus = xPrimeMinus(eta);
  const double _xPrimePlus = xPrimePlus(eta);

  return func(x, _xPrimeMinus, _xPrimePlus, _s, _delta, _B, _psi);
}

NuElectronSpectrum::NuElectronSpectrum()
    : KelnerAharonian2008("data/xsecs_KA2018_nuElectron.txt") {}

double NuElectronSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double NuElectronSpectrum::xPrimeMinus(double eta) const { return xMinus(eta) / 4.; }

double NuElectronSpectrum::xPrimePlus(double eta) const { return xPlus(eta); }

double NuElectronSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xPrimeMinus = xPrimeMinus(eta);
  const double _xPrimePlus = xPrimePlus(eta);

  return func(x, _xPrimeMinus, _xPrimePlus, _s, _delta, _B, _psi);
}

GammaSpectrum::GammaSpectrum() : KelnerAharonian2008("data/xsecs_KA2018_gamma.txt") {}

double GammaSpectrum::psi(double rho) const { return 2.5 + 0.4 * std::log(rho); }

double GammaSpectrum::xPrimeMinus(double eta) const { return 0; }

double GammaSpectrum::xPrimePlus(double eta) const { return 0; }

double GammaSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xMinus = xMinus(eta);
  const double _xPlus = xPlus(eta);

  return func(x, _xMinus, _xPlus, _s, _delta, _B, _psi);
};

ElectronSpectrum::ElectronSpectrum() : KelnerAharonian2008("data/xsecs_KA2018_electron.txt") {}

double ElectronSpectrum::psi(double rho) const {
  return (rho > 4.) ? 6. * (1. - std::exp(1.5 * (4. - rho))) : 0.;
}

double ElectronSpectrum::xPrimeMinus(double eta) const {
  return (eta - 2. * m_r - std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)))) / (4. * (1. + eta));
}

double ElectronSpectrum::xPrimePlus(double eta) const {
  return (eta - 2. * m_r + std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)))) / (2. * (1. + eta));
}

double ElectronSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  if (rho < 2.14) return 0.;

  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xPrimeMinus = xPrimeMinus(eta);
  const double _xPrimePlus = xPrimePlus(eta);

  return func(x, _xPrimeMinus, _xPrimePlus, _s, _delta, _B, _psi);
}

PositronSpectrum::PositronSpectrum() : KelnerAharonian2008("data/xsecs_KA2018_positron.txt") {}

double PositronSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double PositronSpectrum::xPrimeMinus(double eta) const { return 0; }

double PositronSpectrum::xPrimePlus(double eta) const { return 0; }

double PositronSpectrum::Phi(double eta, double x) const {
  if (x > 1.) return 0.;

  const double rho = eta / m_eta_0;  // Eq. 16
  const double _B = B(rho);
  const double _s = s(rho);
  const double _delta = delta(rho);
  const double _psi = psi(rho);
  const double _xMinus = xMinus(eta) / 4.;
  const double _xPlus = xPlus(eta);

  return func(x, _xMinus, _xPlus, _s, _delta, _B, _psi);
}

}  // namespace interactions
