#include "rates.h"

#include "simprop.h"

std::vector<std::vector<double>> loadFromFile(std::string filename, size_t redshiftSize) {
  // Open file
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  // Get file size
  file.seekg(0, std::ios::end);
  std::streamsize fileSize = file.tellg();
  file.seekg(0, std::ios::beg);

  // Check if file size is a multiple of double size
  if (fileSize % sizeof(double) != 0) {
    throw std::runtime_error("Invalid binary file size, not a multiple of double.");
  }

  // Calculate the number of doubles
  size_t numDoubles = fileSize / sizeof(double);

  // Read the data into a vector
  std::vector<double> data(numDoubles);
  if (!file.read(reinterpret_cast<char*>(data.data()), fileSize)) {
    throw std::runtime_error("Error reading binary file.");
  }

  // Close file
  file.close();

  // Slice data in redshift chunks
  std::vector<std::vector<double>> result;
  size_t sliceSize = numDoubles / redshiftSize;
  assert(sliceSize * redshiftSize == numDoubles);
  for (size_t i = 0; i < redshiftSize; i++) {
    result.emplace_back(data.begin() + (i * sliceSize), data.begin() + ((i + 1) * sliceSize - 1));
  }

  return result;
}

namespace gnuprop {

NeutrinoProductionRate::NeutrinoProductionRate() {
  m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
  m_lgx = simprop::utils::LinAxis<double>(m_lgXMin, m_lgXMax, m_xSize);
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
  m_rate = loadFromFile(m_filename, m_zSize);
}

size_t getLowerBoundIndex(double x, const std::vector<double>& v) {
  assert(x >= v.front() and x <= v.back());
  auto it = std::lower_bound(v.begin(), v.end(), x);
  return std::distance(v.begin(), it);
}

double NeutrinoProductionRate::get(double E_nu, double E_p, double z) const {
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

GammaAbsorptionRate::GammaAbsorptionRate() {
  m_z = simprop::utils::LinAxis<double>(m_zMin, m_zMax, m_zSize);
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
  m_rate = loadFromFile(m_filename, m_zSize);
}

double GammaAbsorptionRate::get(double E_gamma, double z) const {
  if (z < m_z.front() || z > m_z.back()) {
    throw std::invalid_argument("z outside the valid range");
  }

  auto iz = getLowerBoundIndex(z, m_z);
  auto array = m_rate[iz];

  const auto lgE = std::log10(E_gamma / SI::eV);

  auto value = 0.;
  if (lgE > m_lgEnergyMin and lgE < m_lgEnergyMax) {
    value += simprop::utils::interpolate(lgE, m_lgE, array);
  }

  return (value * std::pow(1. + z, 3.)) / SI::Gyr;
}

}  // namespace gnuprop