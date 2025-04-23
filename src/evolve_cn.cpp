#include "gnuprop.h"
#include "simprop.h"
#include "tridiag.h"

namespace gnuprop {

double compute_integral(const std::vector<double>& n, const std::vector<double>& E) {
  const auto units = SI::eV / SI::m3;
  double value = 0.;
  for (size_t i = 0; i < n.size(); ++i) {
    value += pow2(E[i]) * n[i];
  }
  return 4. * M_PI / SI::cLight * value / units;
}

void GnuProp::evolve(double zObs) {
  assert(zObs < m_zMax);
  auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
  const auto dz = zAxis[1] - zAxis[0];
  // size_t counter = 0;

  for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {
    auto z = *it;

    LOGD << std::setprecision(5) << z;
    auto IE_p = compute_integral(m_n_proton, m_eAxis);
    auto IE_nu = compute_integral(m_n_nu, m_eAxis);
    auto IE_gamma = compute_integral(m_n_gamma, m_eAxis);
    auto IE_e = compute_integral(m_n_electron, m_eAxis);

    LOGD << "p : " << IE_p;
    LOGD << "n : " << IE_nu;
    LOGD << "g : " << IE_gamma;
    LOGD << "e : " << IE_e;
    LOGD << "test : " << IE_nu / (IE_gamma + IE_e);

    const auto dtdz = std::abs(m_cosmology->dtdz(z));
    const auto H = std::abs(m_cosmology->hubbleRate(z));

    // evolve protons
    {
      evolveProtonEmissivity(z);
      evolveProtonLosses(z);

      std::vector<double> n_proton_up(m_energySize - 1, 0.);
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

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, n_proton_up);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_n_proton[i] = std::max(n_proton_up[i], 0.);
      }
    }

    // evolve neutrinos
    {
      evolveNuEmissivity(z);

      std::vector<double> n_nu_up(m_energySize - 1, 0.);
      for (size_t i = 0; i < m_energySize - 1; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto dE = Eup - E;
        const auto Q = dtdz * m_q_nu[i] - 3. / (1. + z) * m_n_nu[i];
        const auto b = dtdz * E * H;
        const auto bUp = dtdz * Eup * H;
        const auto U_i = 0.5 * dz * bUp / dE;
        const auto C_i = 0.5 * dz * b / dE;
        diagonal[i] = 1. + C_i;
        if (i != m_energySize - 2) upperDiagonal[i] = -U_i;
        knownTerm[i] = U_i * m_n_nu[i + 1] + (1. - C_i) * m_n_nu[i] + dz * Q;
      }

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, n_nu_up);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_n_nu[i] = std::max(n_nu_up[i], 0.);
      }
    }

    // evolve electrons
    {
      evolveElectronEmissivity(z);
      evolveElectronAbsorption(z);

      std::vector<double> n_electron_up(m_energySize - 1, 0.);
      for (size_t i = 0; i < m_energySize - 1; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto dE = Eup - E;
        const auto Q = dtdz * m_q_electron[i] - 3. / (1. + z) * m_n_electron[i];
        const auto b = dtdz * E * H;
        const auto bUp = dtdz * Eup * H;
        const auto U_i = 0.5 * dz * bUp / dE;
        const auto C_i = 0.5 * dz * b / dE;
        diagonal[i] = 1. + C_i;
        if (i != m_energySize - 2) upperDiagonal[i] = -U_i;
        knownTerm[i] = U_i * m_n_electron[i + 1] + (1. - C_i) * m_n_electron[i] + dz * Q;
      }

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, n_electron_up);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        auto C = dtdz * m_k_electron[i];
        auto value = n_electron_up[i] * (1. - 0.5 * dz * C) / (1. + 0.5 * dz * C);
        m_n_electron[i] = std::max(value, 0.);
      }
    }

    // evolve photons
    {
      evolveGammaEmissivity(z);
      evolveGammaAbsorption(z);

      std::vector<double> n_gamma_up(m_energySize - 1, 0.);
      for (size_t i = 0; i < m_energySize - 1; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto dE = Eup - E;
        const auto Q = dtdz * m_q_gamma[i] - 3. / (1. + z) * m_n_gamma[i];
        const auto b = dtdz * E * H;
        const auto bUp = dtdz * Eup * H;
        const auto U_i = 0.5 * dz * bUp / dE;
        const auto C_i = 0.5 * dz * b / dE;
        diagonal[i] = 1. + C_i;
        if (i != m_energySize - 2) upperDiagonal[i] = -U_i;
        knownTerm[i] = U_i * m_n_gamma[i + 1] + (1. - C_i) * m_n_gamma[i] + dz * Q;
      }

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, n_gamma_up);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        auto C = dtdz * m_k_gamma[i];
        auto value = n_gamma_up[i] * (1. - 0.5 * dz * C) / (1. + 0.5 * dz * C);
        m_n_gamma[i] = std::max(value, 0.);
      }
    }
  }
}

}  // namespace gnuprop