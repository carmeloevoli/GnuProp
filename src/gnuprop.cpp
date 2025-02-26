#include "gnuprop.h"

#include <cmath>
#include <fstream>
#include <numeric>

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
  m_losses.push_back(
      std::make_unique<gnuprop::ProtonLossRate>("data/gnuprop_proton_losses_pair.bin"));
  m_losses.push_back(
      std::make_unique<gnuprop::ProtonLossRate>("data/gnuprop_proton_losses_photopion.bin"));

  m_absGammas = std::make_unique<gnuprop::AbsorptionRate>("data/gnuprop_absorption_gammas_cmb.bin");

  m_eAxis = simprop::utils::LogAxis(m_energyMin, m_energyMax, m_energySize);

  m_np.assign(m_energySize, 0.0);
  m_betap.assign(m_energySize, 0.0);

  m_nnu.assign(m_energySize, 0.0);
  m_ngamma.assign(m_energySize, 0.0);

  m_qp.assign(m_energySize, 0.0);
  m_qnu.assign(m_energySize, 0.0);
  m_qgamma.assign(m_energySize, 0.0);

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
        const auto Q = dtdz * std::pow(1. + z, 3.0) * m_qp[i] - 3. / (1. + z) * m_np[i];
        const auto b = dtdz * E * (m_betap[i] + H);
        const auto bUp = dtdz * Eup * (m_betap[i + 1] + H);
        const auto U_i = 0.5 * dz * bUp / dE;
        const auto C_i = 0.5 * dz * b / dE;
        diagonal[i] = 1. + C_i;
        if (i != m_energySize - 2) upperDiagonal[i] = -U_i;
        knownTerm[i] = U_i * m_np[i + 1] + (1. - C_i) * m_np[i] + dz * Q;
      }

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, npUp);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_np[i] = std::max(npUp[i], 0.);
      }
    }

    // evolve neutrinos
    {
      evolveNuEmissivity(z);
      std::vector<double> nnuUp(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf = (m_nnu[i] > 0.) ? (m_nnu[i + 1] / m_nnu[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        nnuUp[i] = m_nnu[i] + dz * (dtdz * m_qnu[i] - m_nnu[i] / (1. + z) * (2. - dlnfdlnE));
      }
      m_nnu = std::move(nnuUp);
    }

    // evolve photons
    {
      evolveGammaEmissivity(z);
      std::vector<double> ngammaUp(m_energySize, 0.);
      // for (size_t i = 0; i < m_energySize - 1.; ++i) {
      //   const auto dlnf = (m_ngamma[i] > 0.) ? (m_ngamma[i + 1] / m_ngamma[i] - 1.) : 0.;
      //   const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
      //   const auto dlnfdlnE = dlnf / dlnE;
      //   const auto K = m_absGammas->get(m_eAxis[i], z);
      //   LOGD << "E: " << m_eAxis[i] << " K / H: " << K / H;
      //   ngammaUp[i] = m_ngamma[i] +
      //                 dz * (dtdz * m_qgamma[i] - m_ngamma[i] / (1. + z) * (2. + K / H -
      //                 dlnfdlnE));
      // }
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf = (m_ngamma[i] > 0.) ? (m_ngamma[i + 1] / m_ngamma[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        const auto K = m_absGammas->get(m_eAxis[i], z);
        const auto C = 1. / (1. + z) * (2. + K / H - dlnfdlnE);
        ngammaUp[i] =
            (dz * dtdz * m_qgamma[i] + (1. - 0.5 * C * dz) * m_ngamma[i]) / (1. + 0.5 * C * dz);
      }

      m_ngamma = std::move(ngammaUp);
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
    m_qp[i] = K * zfactor * std::max(value, 0.0);
  }
}

void GnuProp::evolveProtonLosses(double z) {
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E = m_eAxis[i];
    auto value = 0.;
    for (const auto& loss : m_losses) value += loss->beta(E, z);
    m_betap[i] = value;
  }
}

void GnuProp::evolveNuEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_nu = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_nu : m_photoPionNus) R_j += R_nu->get(E_nu, E_p, z);
      value += m_np[j] * R_j;
    }
    m_qnu[i] = ln_eRatio * value;
  }
}

void GnuProp::evolveGammaEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  for (size_t i = 0; i < m_energySize; ++i) {
    double value = 0.;
    const auto E_gamma = m_eAxis[i];
    for (size_t j = i; j < m_energySize; ++j) {
      const auto E_p = m_eAxis[j];
      double R_j = 0.;
      for (const auto& R_nu : m_photoPionGammas) R_j += R_nu->get(E_gamma, E_p, z);
      value += m_np[j] * R_j;
    }
    m_qgamma[i] = ln_eRatio * value;
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
  const double n2i = SI::cLight / (4.0 * M_PI);

  out << "# E [eV] - I_p - I_nu (1f) - I_gamma [eV-1 m-2 s-1 sr-1] - Q_nu - Q_gamma [eV-1 m-3 "
         "s-1]\n";
  out << std::scientific << std::setprecision(4);

  for (size_t i = 0; i < m_energySize; ++i) {
    auto E = m_eAxis[i] / SI::eV;
    auto I_p = n2i * m_np[i] / I_units;
    auto I_nu = n2i * m_nnu[i] / I_units;
    auto I_gamma = n2i * m_ngamma[i] / I_units;
    auto Q_nu = m_qnu[i] / q_units;
    auto Q_gamma = m_qgamma[i] / q_units;
    out << E << "\t";
    out << I_p << "\t" << I_nu / nFlavors << "\t" << I_gamma << "\t";
    out << Q_nu << "\t" << Q_gamma << "\t";
    out << "\n";
  }

  LOGD << "Dumped spectrum to " << filename;
}

}  // namespace gnuprop
