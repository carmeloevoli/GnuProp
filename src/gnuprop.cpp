#include "gnuprop.h"

#include <omp.h>

#include <cmath>
#include <fstream>

#include "interactions/PhotoPair.h"
#include "simprop.h"

namespace gnuprop {

// Constructor with cosmology only
GnuProp::GnuProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology)
    : m_cosmology(std::move(cosmology)) {
  LOGI << "h: " << m_cosmology->h << " Omega_M: " << m_cosmology->OmegaM
       << " Omega_L: " << m_cosmology->OmegaL;
}

void GnuProp::build() {
  m_sourceNormalization = m_sourceComovingEmissivity * std::abs(m_injSlope - 2.0) / pow2(E_0);
  if (std::fabs(m_injSlope - 2.0) < 1e-3)
    m_sourceNormalization *= 0;  // TO BE DONE
  else if (m_injSlope > 2.)
    m_sourceNormalization *= std::pow(m_leCutoff / E_0, m_injSlope - 2.0);
  else
    m_sourceNormalization *= std::pow(m_heCutoff / E_0, m_injSlope - 2.0);

  m_eAxis = simprop::utils::LogAxis(m_energyMin, m_energyMax, m_energySize);

  m_f_proton.assign(m_energySize, 0.0);
  m_q_proton.assign(m_energySize, 0.0);
  m_beta_proton.assign(m_energySize, 0.0);

  m_f_nu.assign(m_energySize, 0.0);
  m_q_nu.assign(m_energySize, 0.0);

  m_f_gamma.assign(m_energySize, 0.0);
  m_q_gamma.assign(m_energySize, 0.0);
  m_k_gamma.assign(m_energySize, 0.0);

  m_f_electron.assign(m_energySize, 0.0);
  m_q_electron.assign(m_energySize, 0.0);
  m_k_electron.assign(m_energySize, 0.0);

  knownTerm.assign(m_energySize - 1, 0.0);
  diagonal.assign(m_energySize - 1, 0.0);
  upperDiagonal.assign(m_energySize - 2, 0.0);
  lowerDiagonal.assign(m_energySize - 2, 0.0);
}

void GnuProp::computeProtonEmissivity(double z) {
  std::vector<double> q_proton(m_energySize, 0.);
  const double zfactor = (z <= m_zMax) ? std::pow(1.0 + z, m_evolutionIndex) : 0.;

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E = m_eAxis[i];
    auto value = std::pow(E / E_0, -m_injSlope);
    value *= (m_heCutoff > 0.0) ? std::exp(-E / m_heCutoff) : 1.;
    value *= (m_leCutoff > 0.0) ? std::exp(-m_leCutoff / E) : 1.;
    q_proton[i] = m_sourceNormalization * zfactor * std::max(value, 0.0);
  }

  m_q_proton = std::move(q_proton);
}

void GnuProp::computeProtonLosses(double z) {
  std::vector<double> beta_proton(m_energySize, 0.);

#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < m_energySize; ++i) {
    const double E = m_eAxis[i];
    double value = 0.0;
    for (auto const& loss : m_losses_proton) {
      value += loss->beta(E, z);
    }
    beta_proton[i] = value;
  }

  m_beta_proton = std::move(beta_proton);
}

void GnuProp::computeNuEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_nu_pion(m_energySize, 0.);

// Photopion
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_nu = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& q : m_sources_nu_photopion) R_j += q->get(E_nu / E_p, E_p, z);
      value += m_f_proton[j] * R_j;
    }
    q_nu_pion[i] += ln_eRatio * value;
  }

  m_q_nu = std::move(q_nu_pion);
}

void GnuProp::computeGammaEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_g_pion(m_energySize, 0.);

  // Photopion
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_gamma = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& q : m_sources_gamma_photopion) R_j += q->get(E_gamma / E_p, E_p, z);
      value += m_f_proton[j] * R_j;
    }
    q_g_pion[i] += ln_eRatio * value;
  }
  // Inverse Compton
  // #pragma omp parallel for schedule(static)
  //   for (size_t i = 0; i < m_energySize; ++i) {
  //     double value = 0.;
  //     const auto E_gamma = m_eAxis[i];
  //     for (size_t j = i; j < m_energySize; ++j) {
  //       const auto E_electron = m_eAxis[j];
  //       double R_j = 0.;
  //       for (const auto& q : m_sources_gamma_ic) R_j += q->get(E_gamma, E_electron, z);
  //       value += n_electron[j] * R_j;
  //     }
  //     q_gamma[i] += ln_eRatio * value;
  //   }
  m_q_gamma = std::move(q_g_pion);
}

void GnuProp::computeGammaAbsorption(double z) {
  std::vector<double> k_gamma(m_energySize, 0.);

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E_gamma = m_eAxis[i];
    double value = 0.;
    for (const auto& absorption : m_absorption_gamma) {
      value += absorption->get(E_gamma, z);
    }
    k_gamma[i] = value;
  }

  m_k_gamma = std::move(k_gamma);
}

void GnuProp::computeElectronEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_electron(m_energySize, 0.);

// Photopion
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_electron = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& q : m_sources_electron_photopion) R_j += q->get(E_electron / E_p, E_p, z);
      value += m_f_proton[j] * R_j;
    }
    q_electron[i] += ln_eRatio * value;
  }
  // Photopair
  // #pragma omp parallel for schedule(static)
  //   for (size_t i = 0; i < m_energySize; ++i) {
  //     double value = 0.;
  //     const auto E_electron = m_eAxis[i];
  //     for (size_t j = i; j < m_energySize; ++j) {
  //       const auto E_gamma = m_eAxis[j];
  //       double R_j = 0;
  //       for (const auto& q : m_sources_electron_photopair)
  //         R_j += q->get(1. - E_electron / E_gamma, E_gamma, z);
  //       value += m_f_gamma[j] * R_j;
  //     }
  //     q_electron[i] += ln_eRatio * 2.0 * value;
  //   }
  auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 1000);
  const auto ln_xRatio = std::log(xAxis[1] / xAxis[0]);

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_electron = m_eAxis[i];
    for (size_t j = 0; j < xAxis.size(); ++j) {
      const auto x = xAxis[j];
      const auto E_gamma = E_electron / (1. - x);
      double R_j = 0;
      for (const auto& q : m_sources_electron_photopair) R_j += q->get(x, E_gamma, z);
      const auto f_gamma = simprop::utils::interpolate(E_gamma, m_eAxis, m_f_gamma);

      value += x / (1. - x) * R_j * f_gamma;
    }
    q_electron[i] += ln_xRatio * 2.0 * value;
  }

  // const auto E_electron = m_eAxis[i];
  // double value = 0.;
  // for (size_t j = i; j < m_energySize - 1; ++j) {
  //   double R_lo = 0.;
  //   double R_up = 0.;
  //   for (const auto& q : m_sources_electron_photopair) {
  //     R_lo += q->get(1. - E_electron / m_eAxis[j], m_eAxis[j], z);
  //     R_up += q->get(1. - E_electron / m_eAxis[j + 1], m_eAxis[j + 1], z);
  //   }
  //   value += n_gamma[j] * R_lo + n_gamma[j + 1] * R_up;
  // }
  // q_electron[i] += ln_eRatio / 2.0 * 2.0 * value;

  // const auto E_electron = m_eAxis[i];
  // auto integrand = [&](double lnEgamma) {
  //   const auto E_gamma = std::exp(lnEgamma);
  //   double R_j = 0.;
  //   for (const auto& q : m_sources_electron_photopair)
  //     R_j += q->get(1. - E_electron / E_gamma, E_gamma, z);
  //   const auto n_j = simprop::utils::interpolate(E_gamma, m_eAxis, n_gamma);
  //   return n_j * R_j;
  // };
  // const auto a = std::log(E_electron);
  // const auto b = std::log(std::min(E_electron * 1e3, m_energyMax));
  // const size_t N = 10;
  // auto value = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
  // q_electron[i] += 2.0 * value;
  //}
  // Inverse Compton
  // #pragma omp parallel for schedule(static)
  //   for (size_t i = 0; i < m_energySize; ++i) {
  //     double value = 0.;
  //     const auto E_electron = m_eAxis[i];
  //     for (size_t j = i; j < m_energySize; ++j) {
  //       const auto E_electron_prime = m_eAxis[j];
  //       double R_j = 0.;
  //       for (const auto& q : m_sources_electrons_ic) R_j += q->get(E_electron,
  //       E_electron_prime, z); value += n_electron[j] * R_j;
  //     }
  //     q_electron[i] += ln_eRatio * value;
  //   }

  m_q_electron = std::move(q_electron);
}

void GnuProp::computeElectronAbsorption(double z) {
  std::vector<double> k_electron(m_energySize, 0.);

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E_electron = m_eAxis[i];
    double value = 0.;
    for (const auto& absorption : m_absorption_electron) {
      value += absorption->get(E_electron, z);
    }
    k_electron[i] = value;
  }

  m_k_electron = std::move(k_electron);
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
    auto I_p = N2I * m_f_proton[i] / I_units;
    auto I_nu = N2I * m_f_nu[i] / I_units / nFlavors;
    auto I_gamma = N2I * m_f_gamma[i] / I_units;
    auto I_electron = N2I * m_f_electron[i] / I_units;
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
