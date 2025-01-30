#include "losses.h"

#include <fstream>

#include "simprop.h"

std::vector<double> loadFromFile1D(std::string filename) {
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

  file.close();

  return data;
}

namespace gnuprop {

EnergyLosses::EnergyLosses() {
  m_lgE = simprop::utils::LinAxis<double>(m_lgEnergyMin, m_lgEnergyMax, m_energySize);
  m_beta_pair = loadFromFile1D(m_filename_pair);
  m_beta_photopion = loadFromFile1D(m_filename_photopi);
}

double EnergyLosses::beta(double E, double z) const {
  const auto lgE = std::log10(E * (1. + z) / SI::eV);
  auto value = 0.;
  if (lgE > m_lgEnergyMin and lgE < m_lgEnergyMax) {
    value += simprop::utils::interpolate(lgE, m_lgE, m_beta_pair);
    value += (doPhotoPion) ? simprop::utils::interpolate(lgE, m_lgE, m_beta_photopion) : 0;
  }
  return std::pow(1. + z, 3.) * value / SI::Gyr;
}

}  // namespace gnuprop