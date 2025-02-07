#include "rates.h"

#include "cached/load.h"
#include "simprop.h"

namespace gnuprop {

ProtonLossRate::ProtonLossRate(const std::string& filename) {
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
  m_beta = cache::loadFromFile1D(filename);
}

double ProtonLossRate::beta(double E, double z) const {
  const auto lgE = std::log10(E * (1. + z) / SI::eV);
  auto value = 0.;
  if (lgE > m_lgEnergyMin and lgE < m_lgEnergyMax) {
    value = simprop::utils::interpolate(lgE, m_lgE, m_beta);
  }
  return std::pow(1. + z, 3.) * value / SI::Gyr;
}

PhotoPionProductionRate::PhotoPionProductionRate(const std::string& filename) {
  m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
  m_lgx = simprop::utils::LinAxis<double>(m_lgXMin, m_lgXMax, m_xSize);
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
  m_rate = cache::loadFromFile(filename, m_zSize);
}

size_t getLowerBoundIndex(double x, const std::vector<double>& v) {
  assert(x >= v.front() and x <= v.back());
  auto it = std::lower_bound(v.begin(), v.end(), x);
  return std::distance(v.begin(), it);
}

double PhotoPionProductionRate::get(double E_nu, double E_p, double z) const {
  if (z < m_z.front() || z > m_z.back()) {
    throw std::invalid_argument("z outside the valid range");
  }

  auto iz = getLowerBoundIndex(z, m_z);
  auto array = m_rate[iz];

  const auto lgE = std::log10(E_p / SI::eV);
  const auto lgx = std::log10(E_nu / E_p);

  auto value = 0.;
  if (lgE > m_lgEnergyMin and lgE < m_lgEnergyMax and lgx > m_lgXMin and lgx < m_lgXMax) {
    value += simprop::utils::interpolate2d(lgx, lgE, m_lgx, m_lgE, array);
  }

  return (value * std::pow(1. + z, 3.)) / SI::Gyr;
}

// GammaAbsorptionRate::GammaAbsorptionRate() {
//   m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
//   m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
//   m_rate = loadFromFile(m_filename, m_zSize);
// }

// double GammaAbsorptionRate::get(double E_gamma, double z) const {
//   if (z < m_z.front() || z > m_z.back()) {
//     throw std::invalid_argument("z outside the valid range");
//   }

//   auto iz = getLowerBoundIndex(z, m_z);
//   auto array = m_rate[iz];

//   const auto lgE = std::log10(E_gamma / SI::eV);

//   auto value = 0.;
//   if (lgE > m_lgEnergyMin and lgE < m_lgEnergyMax) {
//     value += simprop::utils::interpolate(lgE, m_lgE, array);
//   }

//   return (value * std::pow(1. + z, 3.)) / SI::Gyr;
// }

}  // namespace gnuprop