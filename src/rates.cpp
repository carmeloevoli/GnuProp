#include "rates.h"

#include "cached/load.h"
#include "interpolate.h"
#include "simprop.h"

namespace gnuprop {

// ProtonLossRate
ProtonLossRate::ProtonLossRate(const std::string& filename) {
  auto data = cache::loadFromFile(filename);

  m_lgEnergyMin = std::log10(data[0] / SI::eV);
  m_lgEnergyMax = std::log10(data[1] / SI::eV);
  m_lgEnergySize = static_cast<size_t>(data[2]);

  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_lgEnergySize);
  m_beta = std::vector<double>(data.begin() + 3, data.end());

#ifdef DEBUG
  LOGD << "lgE: " << m_lgEnergyMin << " - " << m_lgEnergyMax << " - " << m_lgEnergySize;
  assert(m_lgEnergyMin < m_lgEnergyMax);
  assert(m_lgEnergyMin > SI::GeV);
  assert(m_lgEnergyMax < 1e16 * SI::GeV);
  assert(m_beta.size() == m_lgEnergySize);
#endif
}

double ProtonLossRate::beta(double E, double z) const {
  const auto lgE = std::log10(E * (1. + z) / SI::eV);

  auto value = 0.;
  if (lgE >= m_lgEnergyMin && lgE <= m_lgEnergyMax) {
    value = utils::interpolate1D(lgE, m_lgEnergyMin, m_lgEnergyMax, m_beta);
  }
  return std::pow(1. + z, 3.) * value / SI::Gyr;
}

// AbsorptionRate
AbsorptionRate::AbsorptionRate(const std::string& filename) {
  auto data = cache::loadFromFile(filename);

  m_zMin = data[0];
  m_zMax = data[1];
  m_zSize = static_cast<size_t>(data[2]);
  m_lgEnergyMin = std::log10(data[3] / SI::eV);
  m_lgEnergyMax = std::log10(data[4] / SI::eV);
  m_lgEnergySize = static_cast<size_t>(data[5]);

  m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_lgEnergySize);

  auto rates = std::vector<double>(data.begin() + 6, data.end());

  size_t sliceSize = m_lgEnergySize;
  assert(sliceSize * m_zSize == rates.size());
  m_rate.reserve(m_zSize);
  for (size_t i = 0; i < m_zSize; i++) {
    m_rate.emplace_back(rates.begin() + (i * sliceSize), rates.begin() + ((i + 1) * sliceSize));
  }

#ifdef DEBUG
  LOGD << "z: " << m_zMin << " - " << m_zMax << " - " << m_zSize;
  LOGD << "lgE: " << m_lgEnergyMin << " - " << m_lgEnergyMax << " - " << m_lgEnergySize;
  assert(m_lgEnergyMin < m_lgEnergyMax);
  assert(m_rate.size() == m_zSize);
  assert(m_rate[0].size() == m_lgEnergySize);
#endif
}

double AbsorptionRate::get(double energy, double z) const {
  if (z < m_zMin || z > m_zMax) {
    throw std::invalid_argument("z outside the valid range");
  }

  const size_t i = utils::getIndex(z, m_zMin, m_zMax, m_zSize);
  const std::vector<double>& array = m_rate[i];
  const auto lgE = std::log10(energy / SI::eV);

  auto value = 0.;
  if (lgE >= m_lgEnergyMin && lgE <= m_lgEnergyMax) {
    value = utils::interpolate1D(lgE, m_lgEnergyMin, m_lgEnergyMax, array);
    value *= std::pow(1. + z, 3.);
  }

  return value / SI::Gyr;
}

// ProductionRate
ProductionRate::ProductionRate(const std::string& filename) {
  auto data = cache::loadFromFile(filename);

  setLimits(data);

  m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
  m_lgx = simprop::utils::LinAxis<double>(m_lgxMin, m_lgxMax, m_xSize);
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_lgEnergySize);

  auto it = data.begin() + 9;

  size_t sliceSize = m_xSize * m_lgEnergySize;
  assert(sliceSize * m_zSize == data.size() - 9);
  m_rate.reserve(m_zSize);
  for (size_t i = 0; i < m_zSize; i++) {
    m_rate.emplace_back(it + (i * sliceSize), it + ((i + 1) * sliceSize));
  }

#ifdef DEBUG
  assert(m_rate.size() == m_zSize);
  assert(m_rate[0].size() == m_xSize * m_lgEnergySize);
#endif
}

void ProductionRate::setLimits(std::vector<double>& data) {
  m_zMin = data[0];
  m_zMax = data[1];
  m_zSize = static_cast<size_t>(data[2]);
  m_lgxMin = std::log10(data[3]);
  m_lgxMax = std::log10(data[4]);
  m_xSize = static_cast<size_t>(data[5]);
  m_lgEnergyMin = std::log10(data[6] / SI::eV);
  m_lgEnergyMax = std::log10(data[7] / SI::eV);
  m_lgEnergySize = static_cast<size_t>(data[8]);
#ifdef DEBUG
  LOGD << "z: " << m_zMin << " - " << m_zMax << " - " << m_zSize;
  LOGD << "lgx: " << m_lgxMin << " - " << m_lgxMax << " - " << m_xSize;
  LOGD << "lgE: " << m_lgEnergyMin << " - " << m_lgEnergyMax << " - " << m_lgEnergySize;
  assert(m_lgEnergyMin < m_lgEnergyMax);
  assert(m_lgxMin < m_lgxMax);
#endif
}

double ProductionRate::get(double E_secondary, double E_primary, double z) const {
  if (z < m_zMin || z > m_zMax) {
    // throw std::invalid_argument("z outside the valid range");
    return 0.;
  }

  size_t i = utils::getIndex(z, m_zMin, m_zMax, m_zSize);
  double lgE = std::log10(E_primary / SI::eV);
  double lgx = std::log10(E_secondary / E_primary);
  const std::vector<double>& array = m_rate[i];

  auto value = 0.;
  if (lgE >= m_lgEnergyMin && lgE <= m_lgEnergyMax && lgx >= m_lgxMin && lgx <= m_lgxMax) {
    value = utils::interpolate2D(lgx, lgE, m_lgxMin, m_lgxMax, m_xSize, m_lgEnergyMin,
                                 m_lgEnergyMax, m_lgEnergySize, array);
    value *= std::pow(1. + z, 3.);
  }

  return value / SI::Gyr;
}

}  // namespace gnuprop