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

    evolveProtonEmissivity(z);
    evolveProtonLosses(z);
    evolveNuEmissivity(z);
    evolveElectronEmissivity(z);
    evolveElectronAbsorption(z);
    evolveGammaEmissivity(z);
    evolveGammaAbsorption(z);

    // evolve protons
    {
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
      std::vector<double> n_nu_up(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf = (m_n_nu[i] > 0.) ? (m_n_nu[i + 1] / m_n_nu[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        const auto C = 1. / (1. + z) * (2. - dlnfdlnE);
        n_nu_up[i] =
            (dz * dtdz * m_q_nu[i] + (1. - 0.5 * C * dz) * m_n_nu[i]) / (1. + 0.5 * C * dz);
      }

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_n_nu[i] = std::max(n_nu_up[i], 0.);
      }
    }

    // evolve electrons
    {
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

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_n_electron[i] = std::max(n_electron_up[i], 0.);
      }
    }

    // evolve photons
    {
      std::vector<double> n_gamma_up(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto dlnf = (m_n_gamma[i] > 0.) ? (m_n_gamma[i + 1] / m_n_gamma[i] - 1.) : 0.;
        const auto dlnE = (m_eAxis[i + 1] / m_eAxis[i] - 1.);
        const auto dlnfdlnE = dlnf / dlnE;
        const auto C = 1. / (1. + z) * (2. + m_k_gamma[i] / H - dlnfdlnE);
        n_gamma_up[i] =
            (dz * dtdz * m_q_gamma[i] + (1. - 0.5 * C * dz) * m_n_gamma[i]) / (1. + 0.5 * C * dz);
      }

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_n_gamma[i] = std::max(n_gamma_up[i], 0.);
      }
    }

    // if (counter % 100 == 0) {
    //   std::string filename = "gnuprop_spectrum_" + std::to_string(counter) + ".txt";
    //   dump(filename);
    // }
    // counter++;
  }
}

}  // namespace gnuprop