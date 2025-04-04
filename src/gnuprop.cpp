#include "gnuprop.h"

#include <omp.h>

#include <cmath>
#include <fstream>

#include "simprop.h"

namespace gnuprop {

// Constructor with cosmology only
GnuProp::GnuProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology)
    : m_cosmology(std::move(cosmology)) {
  LOGI << "h: " << m_cosmology->h << " Omega_M: " << m_cosmology->OmegaM
       << " Omega_L: " << m_cosmology->OmegaL;
}

void GnuProp::build() {
  m_proton_losses.push_back(
      std::make_unique<gnuprop::ProtonLossRate>("tables/gnuprop_proton_losses_pair.bin"));
  m_proton_losses.push_back(
      std::make_unique<gnuprop::ProtonLossRate>("tables/gnuprop_proton_losses_photopion.bin"));
  m_gamma_absorption =
      std::make_unique<gnuprop::AbsorptionRate>("tables/gnuprop_absorption_gammas_cmb.bin");
  m_electron_absorption =
      std::make_unique<gnuprop::AbsorptionRate>("tables/gnuprop_absorption_pairs_cmb.bin");

  m_eAxis = simprop::utils::LogAxis(m_energyMin, m_energyMax, m_energySize);

  m_n_proton.assign(m_energySize, 0.0);
  m_q_proton.assign(m_energySize, 0.0);
  m_beta_proton.assign(m_energySize, 0.0);

  m_n_nu.assign(m_energySize, 0.0);
  m_q_nu.assign(m_energySize, 0.0);

  m_n_gamma.assign(m_energySize, 0.0);
  m_q_gamma.assign(m_energySize, 0.0);
  m_k_gamma.assign(m_energySize, 0.0);

  m_n_electron.assign(m_energySize, 0.0);
  m_q_electron.assign(m_energySize, 0.0);
  m_k_electron.assign(m_energySize, 0.0);

  knownTerm.assign(m_energySize - 1, 0.0);
  diagonal.assign(m_energySize - 1, 0.0);
  upperDiagonal.assign(m_energySize - 2, 0.0);
  lowerDiagonal.assign(m_energySize - 2, 0.0);
}

void GnuProp::evolveProtonEmissivity(double z) {
#ifdef DEBUG
  assert(std::abs(m_injSlope - 2.) > 1e-6);
#endif
  static const double E_0 = 1e18 * SI::eV;
  double K = m_sourceComovingEmissivity * std::abs(m_injSlope - 2.0) / pow2(E_0);
  if (m_injSlope > 2.)
    K *= std::pow(m_leCutoff / E_0, m_injSlope - 2.0);
  else
    K *= std::pow(m_heCutoff / E_0, m_injSlope - 2.0);
  const double zfactor = (z <= m_zMax) ? std::pow(1.0 + z, m_evolutionIndex) : 0.;
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E = m_eAxis[i];
    auto value = std::pow(E / E_0, -m_injSlope);
    value *= (m_heCutoff > 0.0) ? std::exp(-E / m_heCutoff) : 1.;
    value *= (m_leCutoff > 0.0) ? std::exp(-m_leCutoff / E) : 1.;
    m_q_proton[i] = K * zfactor * std::max(value, 0.0);
  }
}

void GnuProp::evolveProtonLosses(double z) {
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E = m_eAxis[i];
    auto value = 0.;
    for (const auto& loss : m_proton_losses) value += loss->beta(E, z);
    m_beta_proton[i] = value;
  }
}

void GnuProp::evolveNuEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_nu(m_energySize, 0.);
// Photopion
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_nu = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_nu : m_nu_photopion) R_j += R_nu->get(E_nu, E_p, z);
      value += m_n_proton[j] * R_j;
    }
    q_nu[i] += ln_eRatio * value;
  }
  m_q_nu = std::move(q_nu);
}

void GnuProp::evolveGammaEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_gamma(m_energySize, 0.);
// Photopion
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_gamma = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_gamma : m_gamma_photopion) R_j += R_gamma->get(E_gamma, E_p, z);
      value += m_n_proton[j] * R_j;
    }
    q_gamma[i] += ln_eRatio * value;
  }
// Inverse Compton
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_gamma = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_electron = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_gamma : m_gamma_ic) R_j += R_gamma->get(E_gamma, E_electron, z);
      value += m_n_electron[j] * R_j;
    }
    q_gamma[i] += ln_eRatio * value;
  }
  m_q_gamma = std::move(q_gamma);
}

void GnuProp::evolveGammaAbsorption(double z) {
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E_gamma = m_eAxis[i];
    m_k_gamma[i] = m_gamma_absorption->get(E_gamma, z);
  }
}

void GnuProp::evolveElectronEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_electron(m_energySize, 0.);
// Photopion
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_electron = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_electron : m_electron_photopion)
        R_j += R_electron->get(E_electron, E_p, z);
      value += m_n_proton[j] * R_j;
    }
    q_electron[i] += ln_eRatio * value;
  }
// Photopair
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_electron = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_gamma = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_electron : m_electron_photopair)
        R_j += R_electron->get(E_electron, E_gamma, z);
      value += m_n_gamma[j] * R_j;
    }
    q_electron[i] += ln_eRatio * value;
  }
// Inverse Compton
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_electron = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_electron_prime = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_electron : m_electrons_ic)
        R_j += R_electron->get(E_electron, E_electron_prime, z);
      value += m_n_electron[j] * R_j;
    }
    q_electron[i] += ln_eRatio * value;
  }
  m_q_electron = std::move(q_electron);
}

void GnuProp::evolveElectronAbsorption(double z) {
#pragma omp parallel for
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E_electron = m_eAxis[i];
    m_k_electron[i] = m_electron_absorption->get(E_electron, z);
  }
}

void GnuProp::dump(const std::string& filename) const {
  std::ofstream out("output/" + filename);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output file: output/" + filename);
  }

  const double I_units = 1.0 / (SI::eV * SI::m2 * SI::sr * SI::sec);
  const double q_units = 1.0 / (SI::eV * SI::m3 * SI::sec);
  const double nFlavors = 3.0;
  const double N2I = SI::cLight / (4.0 * M_PI);

  out << "# E [eV] ";
  out << "- I_p - I_nu (1f) - I_gamma - I_pair [eV-1 m-2 s-1 sr-1] ";
  out << "- Q_nu - Q_gamma - Q_pair [eV-1 m-3 s-1]\n";

  out << std::scientific << std::setprecision(4);

  for (size_t i = 0; i < m_energySize; ++i) {
    auto E = m_eAxis[i] / SI::eV;
    auto I_p = N2I * m_n_proton[i] / I_units;
    auto I_nu = N2I * m_n_nu[i] / I_units / nFlavors;
    auto I_gamma = N2I * m_n_gamma[i] / I_units;
    auto I_electron = N2I * m_n_electron[i] / I_units;
    auto Q_nu = m_q_nu[i] / q_units;
    auto Q_gamma = m_q_gamma[i] / q_units;
    auto Q_electron = m_q_electron[i] / q_units;

    out << E << "\t";
    out << I_p << "\t" << I_nu << "\t" << I_gamma << "\t" << I_electron << "\t";
    out << Q_nu << "\t" << Q_gamma << "\t" << Q_electron << "\t";
    out << "\n";
  }

  LOGD << "Dumped spectrum to " << filename;
}

}  // namespace gnuprop
