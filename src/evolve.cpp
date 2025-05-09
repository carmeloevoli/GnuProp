// #include "ark45.h"
#include "gnuprop.h"
#include "simprop.h"

#define ABSTOL 0
#define RELTOL 1e-3

namespace gnuprop {

int ode_func_protons(double z, const double y[], double f[], void *params) {
  auto *prop = static_cast<GnuProp *>(params);
  const auto eSize = prop->m_energySize;
  const auto dtdz = std::abs(prop->m_cosmology->dtdz(z));
  const auto H = std::abs(prop->m_cosmology->hubbleRate(z));

  for (size_t i = 0; i < eSize - 1; ++i) {
    const auto E = prop->m_eAxis[i];
    const auto Eup = prop->m_eAxis[i + 1];
    const auto dE = Eup - E;
    const auto b = dtdz * E * (prop->m_beta_proton[i] + H);
    const auto bUp = dtdz * Eup * (prop->m_beta_proton[i + 1] + H);

    f[i] = -dtdz * prop->m_q_proton[i];
    f[i] -= (bUp * y[i + 1] - b * y[i]) / dE;
  }
  f[eSize - 1] = 0.0;

  return GSL_SUCCESS;
}

int ode_jac_protons(double z, const double y[], double *dfdy, double dfdt[], void *params) {
  auto *prop = static_cast<GnuProp *>(params);
  const auto eSize = prop->m_energySize;
  const auto dtdz = std::abs(prop->m_cosmology->dtdz(z));
  const auto H = std::abs(prop->m_cosmology->hubbleRate(z));

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, eSize, eSize);
  gsl_matrix *m = &dfdy_mat.matrix;

  for (size_t i = 0; i < eSize - 1; ++i) {
    const auto E = prop->m_eAxis[i];
    const auto Eup = prop->m_eAxis[i + 1];
    const auto dE = Eup - E;
    const auto b = dtdz * E * (prop->m_beta_proton[i] + H);
    const auto bUp = dtdz * Eup * (prop->m_beta_proton[i + 1] + H);

    gsl_matrix_set(m, i, i, b / dE);
    gsl_matrix_set(m, i, i + 1, -bUp / dE);
  }

  for (size_t i = 0; i < eSize; ++i) {
    dfdt[i] = 0.0;
  }

  return GSL_SUCCESS;
}

int ode_func_neutrinos(double z, const double y[], double f[], void *params) {
  auto *prop = static_cast<GnuProp *>(params);
  const auto eSize = prop->m_energySize;
  const auto dtdz = std::abs(prop->m_cosmology->dtdz(z));
  const auto H = std::abs(prop->m_cosmology->hubbleRate(z));

  for (size_t i = 0; i < eSize - 1; ++i) {
    const auto E = prop->m_eAxis[i];
    const auto Eup = prop->m_eAxis[i + 1];
    const auto dE = Eup - E;
    const auto b = dtdz * E * H;
    const auto bUp = dtdz * Eup * H;

    f[i] = -dtdz * prop->m_q_nu[i];
    f[i] -= (bUp * y[i + 1] - b * y[i]) / dE;
  }
  f[eSize - 1] = 0.0;

  return GSL_SUCCESS;
}

int ode_func_gammas(double z, const double y[], double f[], void *params) {
  auto *prop = static_cast<GnuProp *>(params);
  const auto eSize = prop->m_energySize;
  const auto dtdz = std::abs(prop->m_cosmology->dtdz(z));
  const auto H = std::abs(prop->m_cosmology->hubbleRate(z));

  for (size_t i = 0; i < eSize - 1; ++i) {
    const auto E = prop->m_eAxis[i];
    const auto Eup = prop->m_eAxis[i + 1];
    const auto dE = Eup - E;
    const auto b = dtdz * E * H;
    const auto bUp = dtdz * Eup * H;

    f[i] = -dtdz * (prop->m_q_gamma[i] - prop->m_k_gamma[i] * y[i]);
    f[i] -= (bUp * y[i + 1] - b * y[i]) / dE;
  }
  f[eSize - 1] = 0.0;

  return GSL_SUCCESS;
}

int ode_func_electrons(double z, const double y[], double f[], void *params) {
  auto *prop = static_cast<GnuProp *>(params);
  const auto eSize = prop->m_energySize;
  const auto dtdz = std::abs(prop->m_cosmology->dtdz(z));
  const auto H = std::abs(prop->m_cosmology->hubbleRate(z));

  for (size_t i = 0; i < eSize - 1; ++i) {
    const auto E = prop->m_eAxis[i];
    const auto Eup = prop->m_eAxis[i + 1];
    const auto dE = Eup - E;
    const auto b = dtdz * E * H;
    const auto bUp = dtdz * Eup * H;

    f[i] = -dtdz * prop->m_q_electron[i];
    f[i] -= (bUp * y[i + 1] - b * y[i]) / dE;
  }
  f[eSize - 1] = 0.0;

  return GSL_SUCCESS;
}

int ode_jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
  // double k = *((double *)params);
  // (void)t;
  // gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
  // gsl_matrix *m = &dfdy_mat.matrix;
  // gsl_matrix_set(m, 0, 0, 0.0);
  // gsl_matrix_set(m, 0, 1, 1.0);
  // gsl_matrix_set(m, 1, 0, -(1.0 + 2.0 * k * y[0] * y[1]));
  // gsl_matrix_set(m, 1, 1, k * (1.0 * y[0] * y[0]));
  // /* Autonomous. */
  // dfdt[0] = 0.0;
  // dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

void GnuProp::evolve(double zObs) {
  assert(zObs < m_zMax);
  auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
  m_dz = zAxis[1] - zAxis[0];

  const double f_units = SI::eV / SI::m3;
  const double q_units = SI::eV / SI::m3;

  const double hstart = -1e-15;
  const double epsabs = 1e-30;
  const double epsrel = 1e-3;

  std::vector<std::vector<double>> log;

  for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {
    auto z = *it;

    LOGD << std::setprecision(5) << " z : " << z;

    const auto dtdz = std::abs(m_cosmology->dtdz(z));
    const auto H = std::abs(m_cosmology->hubbleRate(z));

    // LOGD << "dtdz : " << dtdz / SI::Gyr << " H^-1 : " << (1. / H) / SI::Gyr;

    std::vector<double> log_z;
    log_z.push_back(z);
    log_z.push_back(H);
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_f_proton));
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_f_nu));
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_f_gamma));
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_f_electron));

    computeProtonEmissivity(z);
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_q_proton));

    computeProtonLosses(z);
    std::vector<double> beta_f(m_energySize, 0.);
    for (std::size_t i = 0; i < m_energySize; ++i) {
      beta_f[i] = m_beta_proton[i] * m_f_proton[i];
    }
    log_z.push_back(compute_energy_integral<double>(m_eAxis, beta_f));

    computeNuEmissivity(z);
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_q_nu));

    computeGammaEmissivity(z);
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_q_gamma));

    computeElectronEmissivity(z);
    log_z.push_back(compute_energy_integral<double>(m_eAxis, m_q_electron));

    computeGammaAbsorption(z);
    std::vector<double> k_g(m_energySize, 0.);
    for (std::size_t i = 0; i < m_energySize; ++i) {
      k_g[i] = m_k_gamma[i] * m_f_gamma[i];
    }
    log_z.push_back(compute_energy_integral<double>(m_eAxis, k_g));

    computeElectronAbsorption(z);
    std::vector<double> k_e(m_energySize, 0.);
    for (std::size_t i = 0; i < m_energySize; ++i) {
      k_e[i] = m_k_electron[i] * m_f_electron[i];
    }
    log_z.push_back(compute_energy_integral<double>(m_eAxis, k_e));

    // evolve protons
    {
      gsl_odeiv2_system sys;
      sys.function = ode_func_protons;
      sys.jacobian = ode_jac_protons;
      sys.dimension = m_energySize;
      sys.params = this;

      const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
      gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, epsabs, epsrel);

      double zStart = z;
      double zEnd = z - m_dz;

      std::vector<double> y = m_f_proton;

      int status = gsl_odeiv2_driver_apply(d, &zStart, zEnd, y.data());

      m_f_proton.assign(y.begin(), y.end());

      if (status != GSL_SUCCESS) {
        std::cerr << "Error: " << gsl_strerror(status) << std::endl;
        throw std::runtime_error("GSL ODE solver failed");
      }

      gsl_odeiv2_driver_free(d);
    }

    // evolve neutrinos
    {
      gsl_odeiv2_system sys;
      sys.function = ode_func_neutrinos;
      sys.jacobian = ode_jac;
      sys.dimension = m_energySize;
      sys.params = this;

      const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
      gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, epsabs, epsrel);

      double zStart = z;
      double zEnd = z - m_dz;

      std::vector<double> y = m_f_nu;

      int status = gsl_odeiv2_driver_apply(d, &zStart, zEnd, y.data());

      m_f_nu.assign(y.begin(), y.end());

      if (status != GSL_SUCCESS) {
        std::cerr << "Error: " << gsl_strerror(status) << std::endl;
        throw std::runtime_error("GSL ODE solver failed");
      }

      gsl_odeiv2_driver_free(d);
    }

    // evolve gammas
    {
      gsl_odeiv2_system sys;
      sys.function = ode_func_gammas;
      sys.jacobian = ode_jac;
      sys.dimension = m_energySize;
      sys.params = this;

      const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
      gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, epsabs, epsrel);

      double zStart = z;
      double zEnd = z - m_dz;

      std::vector<double> y = m_f_gamma;

      int status = gsl_odeiv2_driver_apply(d, &zStart, zEnd, y.data());

      m_f_gamma.assign(y.begin(), y.end());

      if (status != GSL_SUCCESS) {
        std::cerr << "Error: " << gsl_strerror(status) << std::endl;
        throw std::runtime_error("GSL ODE solver failed");
      }

      gsl_odeiv2_driver_free(d);
    }

    // evolve electrons
    {
      gsl_odeiv2_system sys;
      sys.function = ode_func_electrons;
      sys.jacobian = ode_jac;
      sys.dimension = m_energySize;
      sys.params = this;

      const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
      gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, epsabs, epsrel);

      double zStart = z;
      double zEnd = z - m_dz;

      std::vector<double> y = m_f_electron;

      int status = gsl_odeiv2_driver_apply(d, &zStart, zEnd, y.data());

      m_f_electron.assign(y.begin(), y.end());

      if (status != GSL_SUCCESS) {
        std::cerr << "Error: " << gsl_strerror(status) << std::endl;
        throw std::runtime_error("GSL ODE solver failed");
      }

      gsl_odeiv2_driver_free(d);
    }

    // auto IE_p = compute_integral(m_n_proton, m_eAxis);
    // auto IE_nu = compute_integral(m_n_nu, m_eAxis);
    // auto IE_gamma = compute_integral(m_n_gamma, m_eAxis);
    // auto IE_electron = compute_integral(m_n_electron, m_eAxis);

    // LOGD << "p : " << std::setprecision(3) << IE_p;
    // LOGD << "n : " << std::setprecision(3) << IE_nu;
    // LOGD << "g : " << std::setprecision(3) << IE_gamma;
    // LOGD << "e : " << std::setprecision(3) << IE_electron;

    // auto q_p = computeProtonEmissivity(z);
    // auto q_nu = computeNuEmissivity(z, m_n_proton);
    // auto q_gamma = computeGammaEmissivity(z, m_n_proton, m_n_electron);
    // auto q_electron = computeElectronEmissivity(z, m_n_proton, m_n_gamma, m_n_electron);
    // auto K_gamma = computeGammaAbsorption(z);

    // auto QE_p = compute_integral(q_p, m_eAxis);
    // auto QE_nu = compute_integral(q_nu, m_eAxis);
    // auto QE_gamma = compute_integral(q_gamma, m_eAxis);
    // auto QE_electron = compute_integral(q_electron, m_eAxis);

    // std::vector<double> nk_gamma(m_energySize, 0.);
    // for (size_t i = 0; i < m_energySize; ++i) {
    //   nk_gamma[i] = K_gamma[i] * m_n_gamma[i];
    // }

    // LOGD << "Q_p : " << std::setprecision(3) << QE_p;
    // LOGD << "Q_nu : " << std::setprecision(3) << QE_nu;
    // LOGD << "Q_g : " << std::setprecision(3) << QE_gamma;
    // LOGD << "Q_e : " << std::setprecision(3) << QE_electron / QE_gamma;
    // LOGD << "k_g : " << std::setprecision(3) << compute_integral(nk_gamma, m_eAxis) / QE_gamma;

    // if (counter == 30) {
    //   for (size_t i = 0; i < m_energySize; ++i) {
    //     std::cout << m_eAxis[i] / SI::eV << " " << m_n_proton[i] << " " << m_n_nu[i] << " "
    //               << m_n_gamma[i] << " " << m_n_electron[i] << " " << q_gamma[i] << " "
    //               << q_electron[i] << " " << K_gamma[i] << "\n";
    //   }
    //   exit(0);
    // }

    // evolve protons
    {
      // for (size_t i = 0; i < m_energySize; ++i) {
      //   y[i] = m_n_proton[i];
      //   y[i + m_energySize] = m_n_nu[i];
      //   y[i + 2 * m_energySize] = m_n_gamma[i];
      //   y[i + 3 * m_energySize] = m_n_electron[i];
      // }

      // int status = gsl_odeiv2_driver_apply(d, &zStart, zEnd, y);

      // if (status != GSL_SUCCESS) {
      //   std::cerr << "Error: " << gsl_strerror(status) << std::endl;
      //   throw std::runtime_error("GSL ODE solver failed");
      // }

      // for (size_t i = 0; i < m_energySize; ++i) {
      //   m_n_proton[i] = y[i];
      //   m_n_nu[i] = y[i + m_energySize];
      //   m_n_gamma[i] = y[i + 2 * m_energySize];
      //   m_n_electron[i] = y[i + 3 * m_energySize];
      // }

      //      gsl_odeiv2_driver* d =
      //        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_odeiv2_stepper, 1e-6, 1e-6,
      //        0.1);
      // const auto ode_2d = [&](double z, const std::vector<double>& state) {
      //   evolveProtonEmissivity(z);
      //   evolveProtonLosses(z);
      //   const auto dtdz = std::abs(m_cosmology->dtdz(z));
      //   const auto H = std::abs(m_cosmology->hubbleRate(z));

      //   std::vector<double> dydt(state.size(), 0.);
      //   for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
      //     const auto E = m_eAxis[i];
      //     const auto Eup = m_eAxis[i + 1];
      //     const auto dE = Eup - E;
      //     const auto b = dtdz * E * (m_beta_proton[i] + H);
      //     const auto bUp = dtdz * Eup * (m_beta_proton[i + 1] + H);

      //     dydt.at(i) = -dtdz * m_q_proton[i];
      //     dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
      //   }
      //   return dydt;
      // };
      // AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-12, 0.1);

      // std::vector<double> current_state = m_n_proton;
      // auto result = solver_2d.solve(z, current_state, z - dz);
      // m_n_proton = std::move(result);
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
      // const auto ode_2d = [&](double z, const std::vector<double>& state) {
      //   evolveElectronEmissivity(z);
      //   evolveElectronAbsorption(z);
      //   const auto dtdz = std::abs(m_cosmology->dtdz(z));
      //   const auto H = std::abs(m_cosmology->hubbleRate(z));

      //   std::vector<double> dydt(state.size(), 0.);
      //   for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
      //     const auto E = m_eAxis[i];
      //     const auto Eup = m_eAxis[i + 1];
      //     const auto dE = Eup - E;
      //     const auto b = dtdz * E * H;
      //     const auto bUp = dtdz * Eup * H;

      //     dydt.at(i) = -dtdz * (m_q_electron[i] - m_k_electron[i] * state[i]);
      //     dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
      //   }
      //   return dydt;
      // };
      // AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-20, 0.1);

      // std::vector<double> current_state = m_n_electron;
      // auto result = solver_2d.solve(z, current_state, z - dz);
      // m_n_electron = std::move(result);
    }

    // evolve photons
    {
      // const auto ode_2d = [&](double z, const std::vector<double>& state) {
      //   evolveGammaEmissivity(z);
      //   evolveGammaAbsorption(z);
      //   const auto dtdz = std::abs(m_cosmology->dtdz(z));
      //   const auto H = std::abs(m_cosmology->hubbleRate(z));

      //   std::vector<double> dydt(state.size(), 0.);
      //   for (size_t i = 0; i < m_eAxis.size() - 1; ++i) {
      //     const auto E = m_eAxis[i];
      //     const auto Eup = m_eAxis[i + 1];
      //     const auto dE = Eup - E;
      //     const auto b = dtdz * E * H;
      //     const auto bUp = dtdz * Eup * H;

      //     dydt.at(i) = -dtdz * (m_q_gamma[i] - m_k_gamma[i] * state[i]);
      //     dydt.at(i) -= (bUp * state[i + 1] - b * state[i]) / dE;
      //   }
      //   return dydt;
      // };
      // AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, ABSTOL, RELTOL, 1e-20, 0.1);

      // std::vector<double> current_state = m_n_gamma;
      // auto result = solver_2d.solve(z, current_state, z - dz);
      // m_n_gamma = std::move(result);
    }
  }
}

}  // namespace gnuprop