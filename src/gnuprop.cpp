#include "gnuprop.h"

#include <cmath>
#include <fstream>

#include "simprop.h"
#include "tridiag.h"

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

void GnuProp::evolve(double zObs) {
  assert(zObs < m_zMax);
  auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
  const auto dz = zAxis[1] - zAxis[0];

  for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {
    auto z = *it;
    LOGD << std::setprecision(5) << z;
    LOGD << *std::min_element(m_n_gamma.begin(), m_n_gamma.end()) << " "
         << *std::max_element(m_n_gamma.begin(), m_n_gamma.end());

    const auto dtdz = std::abs(m_cosmology->dtdz(z));
    const auto H = std::abs(m_cosmology->hubbleRate(z));

    // evolve protons
    {
      evolveProtonEmissivity(z);
      evolveProtonLosses(z);
      std::vector<double> npUp(m_energySize - 1, 0.);
      for (size_t i = 0; i < m_energySize - 1; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto dE = Eup - E;
        const auto Q = dtdz * std::pow(1. + z, 3.0) * m_q_proton[i] - 3. / (1. + z) * m_n_proton[i];
        const auto b = dtdz * E * (m_beta_proton[i] + H);
        const auto bUp = dtdz * Eup * (m_beta_proton[i + 1] + H);
        const auto U_i = 0.5 * dz * bUp / dE;
        const auto C_i = 0.5 * dz * b / dE;
        diagonal[i] = 1. + C_i;
        if (i != m_energySize - 2) upperDiagonal[i] = -U_i;
        knownTerm[i] = U_i * m_n_proton[i + 1] + (1. - C_i) * m_n_proton[i] + dz * Q;
      }

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, npUp);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_n_proton[i] = std::max(npUp[i], 0.);
      }
    }

    // evolve neutrinos
    {
      evolveNuEmissivity(z);
      std::vector<double> n_nu_up(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf = (m_n_nu[i] > 0.) ? (m_n_nu[i + 1] / m_n_nu[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        const auto C = 1. / (1. + z) * (2. - dlnfdlnE);
        n_nu_up[i] =
            (dz * dtdz * m_q_nu[i] + (1. - 0.5 * C * dz) * m_n_nu[i]) / (1. + 0.5 * C * dz);
      }
      m_n_nu = std::move(n_nu_up);
    }

    // evolve photons
    {
      evolveGammaEmissivity(z);
      evolveGammaAbsorption(z);
      std::vector<double> n_gamma_up(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf = (m_n_gamma[i] > 0.) ? (m_n_gamma[i + 1] / m_n_gamma[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        const auto C = 1. / (1. + z) * (2. + m_k_gamma[i] / H - dlnfdlnE);
        n_gamma_up[i] =
            (dz * dtdz * m_q_gamma[i] + (1. - 0.5 * C * dz) * m_n_gamma[i]) / (1. + 0.5 * C * dz);
      }
      m_n_gamma = std::move(n_gamma_up);
    }

    // evolve electrons
    {
      evolveElectronEmissivity(z);
      evolveElectronAbsorption(z);
      std::vector<double> n_electron_up(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf =
            (m_n_electron[i] > 0.) ? (m_n_electron[i + 1] / m_n_electron[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        const auto C = 1. / (1. + z) * (2. + m_k_electron[i] / H - dlnfdlnE);
        n_electron_up[i] = (dz * dtdz * m_q_electron[i] + (1. - 0.5 * C * dz) * m_n_electron[i]) /
                           (1. + 0.5 * C * dz);
      }
      m_n_electron = std::move(n_electron_up);
    }
  }
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
  // photopion
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
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_gamma = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_electron = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_gamma : m_gamma_ic) R_j += R_gamma->get(E_gamma, E_electron, z);
      value += m_n_electron[j] * R_j;
    }
    // q_gamma[i] += ln_eRatio * value;
  }
  m_q_gamma = std::move(q_gamma);
}

void GnuProp::evolveGammaAbsorption(double z) {
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E_gamma = m_eAxis[i];
    m_k_gamma[i] = m_gamma_absorption->get(E_gamma, z);
  }
}

void GnuProp::evolveElectronEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> q_electron(m_energySize, 0.);
  // Photopion
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
