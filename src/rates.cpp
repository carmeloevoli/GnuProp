#include "rates.h"

#include "cached/load.h"
#include "interpolate.h"
#include "simprop.h"

namespace gnuprop {

// ProtonLossRate
ProtonLossRate::ProtonLossRate(const std::string& filename) {
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_lgEnergySize);
  m_beta = cache::loadFromFile1D(filename);
#ifdef DEBUG
  assert(m_beta.size() == m_lgEnergySize);
#endif
}

double ProtonLossRate::beta(double E, double z) const {
  const auto lgE = std::log10(E * (1. + z) / SI::eV);
  auto value = 0.;
  if (lgE > m_lgEnergyMin && lgE < m_lgEnergyMax) {
    value = utils::interpolate1D(lgE, m_lgEnergyMin, m_lgEnergyMax, m_beta);
  }
  return std::pow(1. + z, 3.) * value / SI::Gyr;
}

// AbsorptionRate
AbsorptionRate::AbsorptionRate(const std::string& filename) {
  m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_lgEnergySize);
  m_rate = cache::loadFromFile(filename, m_zSize);
#ifdef DEBUG
  assert(m_rate.size() == m_zSize);
  assert(m_rate[0].size() == m_lgEnergySize);
#endif
}

double AbsorptionRate::get(double energy, double z) const {
  if (z < m_zMin || z > m_zMax) {
    throw std::invalid_argument("z outside the valid range");
  }

  const size_t i = utils::getIndex(z, m_zMin, m_zMax, m_zSize);
  // const std::vector<double>& array = m_rate[i];

  const auto lgE = std::log10(energy / SI::eV);

  auto value = 0.;
  if (lgE >= m_lgEnergyMin && lgE <= m_lgEnergyMax) {
    value = utils::interpolate1D(lgE, m_lgEnergyMin, m_lgEnergyMax, m_rate[i]);
    value *= std::pow(1. + z, 3.);
  }

  return value / SI::Gyr;
}

// // PhotoPionRate
// PhotoPionRate::PhotoPionRate(const std::string& filename) {
//   m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
//   m_lgx = simprop::utils::LinAxis<double>(m_lgXMin, m_lgXMax, m_xSize);
//   m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
//   m_rate = cache::loadFromFile(filename, m_zSize);
// #ifdef DEBUG
//   assert(m_rate.size() == m_zSize);
//   assert(m_rate[0].size() == m_xSize * m_energySize);
// #endif
// }

// double PhotoPionRate::get(double E_nu, double E_p, double z) const {
//   if (z < m_zMin || z > m_zMax) {
//     throw std::invalid_argument("z outside the valid range");
//   }

//   const size_t i = utils::getIndex(z, m_zMin, m_zMax, m_zSize);
//   const std::vector<double>& array = m_rate[i];

//   const auto lgE = std::log10(E_p / SI::eV);
//   const auto lgx = std::log10(E_nu / E_p);

//   auto value = 0.;
//   if (lgE > m_lgEnergyMin && lgE < m_lgEnergyMax && lgx > m_lgXMin && lgx < m_lgXMax) {
//     value = utils::interpolate2D(lgx, lgE, m_lgXMin, m_lgXMax, m_xSize, m_lgEnergyMin,
//                                  m_lgEnergyMax, m_energySize, array);
//     value *= std::pow(1. + z, 3.);
//   }

//   return value / SI::Gyr;
// }

// // PhotoPairRate
// PhotoPairRate::PhotoPairRate(const std::string& filename) {
//   m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
//   m_lgx = simprop::utils::LinAxis<double>(m_lgXMin, m_lgXMax, m_xSize);
//   m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
//   m_rate = cache::loadFromFile(filename, m_zSize);
// #ifdef DEBUG
//   assert(m_rate.size() == m_zSize);
//   assert(m_rate[0].size() == m_xSize * m_energySize);
// #endif
// }

// double PhotoPairRate::get(double E_nu, double E_p, double z) const {
//   if (z < m_zMin || z > m_zMax) {
//     throw std::invalid_argument("z outside the valid range");
//   }

//   const size_t i = utils::getIndex(z, m_zMin, m_zMax, m_zSize);
//   const std::vector<double>& array = m_rate[i];

//   const auto lgE = std::log10(E_p / SI::eV);
//   const auto lgx = std::log10(E_nu / E_p);

//   auto value = 0.;
//   if (lgE > m_lgEnergyMin && lgE < m_lgEnergyMax && lgx > m_lgXMin && lgx < m_lgXMax) {
//     value = utils::interpolate2D(lgx, lgE, m_lgXMin, m_lgXMax, m_xSize, m_lgEnergyMin,
//                                  m_lgEnergyMax, m_energySize, array);
//     value *= std::pow(1. + z, 3.);
//   }

//   return value / SI::Gyr;
// }

}  // namespace gnuprop