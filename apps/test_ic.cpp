#include <fstream>

#include "gnuprop.h"
#include "interactions/InverseCompton.h"
#include "simprop.h"

#define LIMIT 1000

#define REDSHIFT 4.5

const auto m_ic = std::make_unique<Interactions::InverseCompton>();
const auto m_phField = std::make_unique<simprop::photonfields::CMB>();
const auto precision = 1e-3;
const auto q_gamma =
    std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_inversecompton_gammas_cmb.bin");

#undef GAUSSIAN

#ifdef GAUSSIAN
double f_electron(double E_gamma) {
  const double N_0 = 1e20;
  const double EeV = 1e18 * SI::eV;
  return N_0 * std::pow(E_gamma / SI::PeV, -2.0) * std::exp(-E_gamma / EeV);
}
#else
double f_electron(double E_electron) {
  const double E_0 = 1e18 * SI::eV;
  const double sigma = 0.2 * E_0;
  return 1e10 * std::exp(-pow2((E_electron - E_0) / sigma));
}
#endif

double K(double E_electron, double z = 0.) {
  auto integrandOuter = [&](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);

    auto integrandInner = [&](double mu) {
      const auto beta = 1.;  // std::sqrt(1. - 4. * SI::electronMassC2 / mu);
      return 0.5 * (1. - beta * mu) * m_ic->sigma_lab(E_electron, epsilon, mu);
    };

    auto value = simprop::utils::QAGIntegration<double>(integrandInner, -1, 1, LIMIT, 1e-5);

    return epsilon * m_phField->density(epsilon, z) * value;
  };

  const auto epsThr = 0.;  // TBD
  const auto epsMin = std::max(epsThr, m_phField->getMinPhotonEnergy());
  const auto epsMax = m_phField->getMaxPhotonEnergy();

  auto value = simprop::utils::QAGIntegration<double>(integrandOuter, std::log(epsMin),
                                                      std::log(epsMax), LIMIT, precision);
  value *= SI::cLight;

  return std::max(value, 0.);
}

double K_gnuprop(double E_electron, double z = 0.) {
  const auto me_2 = pow2(SI::electronMassC2);
  const auto beta_e = std::sqrt(1. - me_2 / pow2(E_electron));
  const auto epsMax = m_phField->getMaxPhotonEnergy();

  auto integrandOuter = [&](double s) {
    const auto epsThr = (s - me_2) / 2. / E_electron / (1. + beta_e);
    const auto epsMin = std::max(epsThr, m_phField->getMinPhotonEnergy());

    auto integrandInner = [&](double lnEpsilon) {
      const auto epsilon = std::exp(lnEpsilon);
      return m_phField->density(epsilon, z) / epsilon;
    };

    const auto a = std::log(epsMin);
    const auto b = std::log(epsMax);
    const size_t N = 10000;
    auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, N, 1e-3);

    return (s - me_2) * m_ic->sigma_com(s) * value;
  };

  const auto sMin = me_2;
  const auto sMax = me_2 + 2. * E_electron * epsMax * (1. + beta_e);

  const size_t N = 10000;
  auto value = simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, N, precision);
  value *= SI::cLight / 8. / beta_e / pow2(E_electron);

  return std::max(value, 0.);
}

double R(double E_gamma, double E_electron, double z = 0.) {
  if (E_electron < E_gamma) return 0.;

  const auto me_2 = pow2(SI::electronMassC2);
  const auto epsTh = 0.;  // me_2 / E_gamma;
  const auto epsMin = std::max(epsTh, m_phField->getMinPhotonEnergy());
  const auto epsMax = m_phField->getMaxPhotonEnergy();

  if (epsMin > epsMax) return 0.;

  auto integrand = [&](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    auto value = m_phField->density(epsilon, z) * m_ic->dsigma_dE(E_electron, epsilon, E_gamma);
    return epsilon * value;
  };

  const auto a = std::log(epsMin);
  const auto b = std::log(epsMax);
  const size_t N = 10000;
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, precision);
  value *= SI::cLight * E_electron;

  return std::max(value, 0.);
}

double Q(double E_gamma, double z = 0.) {
  auto integrand = [&](double lnEelectron) {
    const auto E_electron = std::exp(lnEelectron);
    return R(E_gamma, E_electron, z) * f_electron(E_electron);
  };

  const auto a = std::log(E_gamma);
  const auto b = std::log(1e22 * SI::eV);
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, LIMIT, precision);

  return std::max(value, 0.);
}

double Q_interpolated(double E_gamma, double z = 0.) {
  auto integrand = [&](double lnEelectron) {
    const auto E_electron = std::exp(lnEelectron);
    return q_gamma->get(1. - E_gamma / E_electron, E_electron, z) * f_electron(E_electron);
  };

  const auto a = std::log(E_gamma);
  const auto b = std::log(1e22 * SI::eV);
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, LIMIT, 1e-2);

  return std::max(value, 0.);
}

void print_absorption_term(std::string filename = "output/test_electron_absorption.txt") {
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cerr << "Error: could not open output.txt for writing\n";
  }

  ofs << std::scientific << std::setprecision(3);
  for (double E_e = 1e12 * SI::eV; E_e < 1e23 * SI::eV; E_e *= 1.02) {
    ofs << E_e / SI::eV << "   ";
    ofs << pow2(E_e / SI::eV) * K(E_e, REDSHIFT) * f_electron(E_e) << "   ";
    ofs << pow2(E_e / SI::eV) * K_gnuprop(E_e, REDSHIFT) * f_electron(E_e) << "\n";
  }

  ofs.close();

  std::cout << "written to " << filename << "\n";

  auto I_K = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_electron = std::exp(lnE);
        return pow2(E_electron / SI::eV) * K(E_electron, REDSHIFT) * f_electron(E_electron);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);

  std::cout << std::setprecision(3) << "I_K = " << I_K << "\n";
}

void print_source_terms(std::string filename = "output/test_gamma_source.txt") {
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cerr << "Error: could not open output.txt for writing\n";
  }

  ofs << std::scientific << std::setprecision(3);
  for (double E_gamma = 1e12 * SI::eV; E_gamma < 1e23 * SI::eV; E_gamma *= 1.02) {
    ofs << E_gamma / SI::eV << "   ";
    ofs << pow2(E_gamma / SI::eV) * Q(E_gamma, REDSHIFT) << "  ";
    ofs << pow2(E_gamma / SI::eV) * Q_interpolated(E_gamma, REDSHIFT) << "\n";
  }

  ofs.close();

  auto I_Q = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_gamma = std::exp(lnE);
        return pow2(E_gamma / SI::eV) * Q(E_gamma, REDSHIFT);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);
  std::cout << std::setprecision(3) << "I_Q = " << I_Q << "\n";

  auto I_Q_tilde = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_gamma = std::exp(lnE);
        return pow2(E_gamma / SI::eV) * Q_interpolated(E_gamma, REDSHIFT);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);
  std::cout << std::setprecision(3) << "I_Q_tilde = " << I_Q_tilde << "\n";

  std::cout << "written to " << filename << "\n";
}

void test_interpolation(std::string filename = "output/test_ic_interpolation.txt") {
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cerr << "Error: could not open output.txt for writing\n";
  }

  ofs << std::scientific << std::setprecision(3);
  const double E_electron = 1e18 * SI::eV;
  for (double E_gamma = 1e12 * SI::eV; E_gamma < 1e21 * SI::eV; E_gamma *= 1.02) {
    ofs << E_gamma / SI::eV << "   ";
    ofs << R(E_gamma, E_electron) << "   ";
    ofs << q_gamma->get(1. - E_gamma / E_electron, E_electron) << "\n";
  }
  ofs.close();
}

int main() {
  try {
    print_absorption_term();
    print_source_terms();
    test_interpolation();
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what();
  }
  return 0;
}