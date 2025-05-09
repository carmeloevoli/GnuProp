#include "ark45.h"
#include "gnuprop.h"
#include "simprop.h"

#define ABSTOL 0
#define RELTOL 1e-3

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

  for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {
    auto z = *it;

    LOGD << std::setprecision(5) << z;

    auto IE_p = compute_integral(m_n_proton, m_eAxis);
    auto IE_nu = compute_integral(m_n_nu, m_eAxis);
    auto IE_gamma = compute_integral(m_n_gamma, m_eAxis);
    auto IE_e = compute_integral(m_n_electron, m_eAxis);

    LOGD << "p : " << std::setprecision(2) << IE_p;
    LOGD << "n : " << std::setprecision(2) << IE_nu / IE_p;
    LOGD << "g : " << std::setprecision(2) << IE_gamma / IE_p;
    LOGD << "e : " << std::setprecision(2) << IE_e / IE_p;

    //    LOGD << "nu / gamma+pairs : " << IE_nu / (IE_gamma + IE_e);

    // evolve protons
    {
      const auto ode_2d = [&](double z, const std::vector<double>& state) {
        evolveProtonEmissivity(z);
        evolveProtonLosses(z);
        const auto dtdz = std::abs(m_cosmology->dtdz(z));
        const auto H = std::abs(m_cosmology->hubbleRate(z));

        std::vector<double> dydt(state.size(), 0.);
        for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
          const auto E = m_eAxis[i];
          const auto Eup = m_eAxis[i + 1];
          const auto dE = Eup - E;
          const auto b = dtdz * E * (m_beta_proton[i] + H);
          const auto bUp = dtdz * Eup * (m_beta_proton[i + 1] + H);

          dydt.at(i) = -dtdz * m_q_proton[i];
          dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
        }
        return dydt;
      };
      AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-12, 0.1);

      std::vector<double> current_state = m_n_proton;
      auto result = solver_2d.solve(z, current_state, z - dz);
      m_n_proton = std::move(result);
    }

    // evolve neutrinos
    {
      // const auto ode_2d = [&](double z, const std::vector<double>& state) {
      //   evolveNuEmissivity(z);
      //   const auto dtdz = std::abs(m_cosmology->dtdz(z));
      //   const auto H = std::abs(m_cosmology->hubbleRate(z));

      //   std::vector<double> dydt(state.size(), 0.);
      //   for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
      //     const auto E = m_eAxis[i];
      //     const auto Eup = m_eAxis[i + 1];
      //     const auto dE = Eup - E;
      //     const auto b = dtdz * E * H;
      //     const auto bUp = dtdz * Eup * H;

      //     dydt.at(i) = -dtdz * m_q_nu[i];
      //     dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
      //   }
      //   return dydt;
      // };
      // AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-12, 0.1);

      // std::vector<double> current_state = m_n_nu;
      // auto result = solver_2d.solve(z, current_state, z - dz);
      // m_n_nu = std::move(result);
    }

    // evolve electrons
    {
      const auto ode_2d = [&](double z, const std::vector<double>& state) {
        evolveElectronEmissivity(z);
        evolveElectronAbsorption(z);
        const auto dtdz = std::abs(m_cosmology->dtdz(z));
        const auto H = std::abs(m_cosmology->hubbleRate(z));

        std::vector<double> dydt(state.size(), 0.);
        for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
          const auto E = m_eAxis[i];
          const auto Eup = m_eAxis[i + 1];
          const auto dE = Eup - E;
          const auto b = dtdz * E * H;
          const auto bUp = dtdz * Eup * H;

          dydt.at(i) = -dtdz * (m_q_electron[i] - m_k_electron[i] * state[i]);
          dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
        }
        return dydt;
      };
      AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-20, 0.1);

      std::vector<double> current_state = m_n_electron;
      auto result = solver_2d.solve(z, current_state, z - dz);
      m_n_electron = std::move(result);
    }

    // evolve photons
    {
      const auto ode_2d = [&](double z, const std::vector<double>& state) {
        evolveGammaEmissivity(z);
        evolveGammaAbsorption(z);
        const auto dtdz = std::abs(m_cosmology->dtdz(z));
        const auto H = std::abs(m_cosmology->hubbleRate(z));

        std::vector<double> dydt(state.size(), 0.);
        for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
          const auto E = m_eAxis[i];
          const auto Eup = m_eAxis[i + 1];
          const auto dE = Eup - E;
          const auto b = dtdz * E * H;
          const auto bUp = dtdz * Eup * H;

          dydt.at(i) = -dtdz * (m_q_gamma[i] - m_k_gamma[i] * state[i]);
          dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
        }
        return dydt;
      };
      AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-20, 0.1);

      std::vector<double> current_state = m_n_gamma;
      auto result = solver_2d.solve(z, current_state, z - dz);
      m_n_gamma = std::move(result);
    }
  }
}

}  // namespace gnuprop