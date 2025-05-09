#include <fstream>

#include "gnuprop.h"
#include "interactions/PhotoPair.h"
#include "simprop.h"

#define LIMIT 10000
#define REDSHIFT 4.5

const auto m_photoPair = std::make_unique<Interactions::PhotoPair>();
const auto m_phField = std::make_unique<simprop::photonfields::CMB>();
const auto precision = 1e-3;
const auto q_electron =
    std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopair_pairs_cmb.bin");

#undef GAUSSIAN

class F_GAMMA {
 public:
  F_GAMMA() {
    _energies = simprop::utils::LogAxis<double>(1e12 * SI::eV, 1e23 * SI::eV, 1000);
    for (size_t i = 0; i < _energies.size(); ++i) {
      _fluxes.push_back(_f_gamma(_energies[i]));
    }
  };
  ~F_GAMMA() = default;

  double operator()(double E_gamma) const {
    auto it = std::lower_bound(_energies.begin(), _energies.end(), E_gamma);
    if (it == _energies.end()) {
      return _fluxes.back();
    }
    if (it == _energies.begin()) {
      return _fluxes.front();
    }
    size_t index = std::distance(_energies.begin(), it);
    double x1 = std::log(_energies[index - 1]);
    double x2 = std::log(_energies[index]);
    double y1 = _fluxes[index - 1];
    double y2 = _fluxes[index];
    double interpolatedFlux = y1 + (y2 - y1) * (std::log(E_gamma) - x1) / (x2 - x1);
    return interpolatedFlux;
  }

 private:
  std::vector<double> _energies;
  std::vector<double> _fluxes;

#ifdef GAUSSIAN
  double f_gamma(double E_gamma) {
    const double E_0 = 1e16 * SI::eV;
    const double sigma = 0.2 * E_0;
    return 1e20 * std::exp(-pow2((E_gamma - E_0) / sigma));
  }
#else
  double _f_gamma(double E_gamma) {
    const double N_0 = 1e20;
    const double EeV = 1e19 * SI::eV;
    return N_0 * std::pow(E_gamma / SI::PeV, -1.) * std::exp(-E_gamma / EeV);
  }
#endif
};

auto f_gamma = F_GAMMA();

double K(double E_gamma, double z = 0.) {
  auto integrandOuter = [&](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);

    auto integrandInner = [&](double mu) {
      const auto beta = 1.;
      return 0.5 * (1. - beta * mu) * m_photoPair->sigma_lab(E_gamma, epsilon, mu);
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

double R(double E_electron, double E_gamma, double z = 0.) {
  if (E_electron > E_gamma) return 0.;

  const auto me_2 = pow2(SI::electronMassC2);
  const auto epsTh = me_2 / E_gamma;
  const auto epsMin = std::max(epsTh, m_phField->getMinPhotonEnergy());
  const auto epsMax = m_phField->getMaxPhotonEnergy();

  if (epsMin > epsMax) return 0.;

  auto integrand = [&](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    auto value =
        m_phField->density(epsilon, z) * m_photoPair->dsigma_dE(E_gamma, epsilon, E_electron);
    return epsilon * value;
  };

  const auto a = std::log(epsMin);
  const auto b = std::log(epsMax);
  const size_t N = 10000;
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, precision);
  value *= SI::cLight * E_gamma;

  return std::max(value, 0.);
}

double Q(double E_electron, double z = 0.) {
  auto integrand = [&](double lnEGamma) {
    const auto E_gamma = std::exp(lnEGamma);
    return R(E_electron, E_gamma, z) * f_gamma(E_gamma);
  };

  const auto a = std::log(E_electron);
  const auto b = std::log(1e22 * SI::eV);
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, LIMIT, precision);

  return std::max(value, 0.);
}

double Q_interpolated(double E_electron, double z = 0.) {
  auto integrand = [&](double lnEGamma) {
    const auto E_gamma = std::exp(lnEGamma);
    return q_electron->get(1. - E_electron / E_gamma, E_gamma, z) * f_gamma(E_gamma);
  };

  const auto a = std::log(E_electron);
  const auto b = std::log(1e22 * SI::eV);
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, LIMIT, 1e-2);

  return std::max(value, 0.);
}

double Q_interpolated_x(double E_electron, double z = 0.) {
  auto x_vec = simprop::utils::LogAxis<double>(1e-5, 1., 1000);
  auto value = 0.;
  for (size_t i = 0; i < x_vec.size() - 1; ++i) {
    const auto x = x_vec[i];
    const auto x_next = x_vec[i + 1];
    auto dx = x_next - x;
    const auto E_gamma = E_electron / (1. - x);
    value += dx / (1. - x) * q_electron->get(x, E_gamma, z) * f_gamma(E_gamma);
  }

  return std::max(value, 0.);
}

void print_absorption_term(std::string filename = "output/test_gamma_absorption.txt") {
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cerr << "Error: could not open output.txt for writing\n";
  }

  ofs << std::scientific << std::setprecision(3);
  for (double E_gamma = 1e12 * SI::eV; E_gamma < 1e23 * SI::eV; E_gamma *= 1.02) {
    ofs << E_gamma / SI::eV << "   ";
    ofs << pow2(E_gamma / SI::eV) * K(E_gamma, REDSHIFT) * f_gamma(E_gamma) << "   ";
    ofs << K(E_gamma, REDSHIFT) << "\n";
  }

  ofs.close();

  std::cout << "written to " << filename << "\n";

  auto I_K = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_gamma = std::exp(lnE);
        return pow2(E_gamma / SI::eV) * K(E_gamma, REDSHIFT) * f_gamma(E_gamma);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);

  std::cout << std::setprecision(3) << "I_K = " << I_K << "\n";
}

void print_source_terms(std::string filename = "output/test_electron_source.txt") {
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cerr << "Error: could not open output.txt for writing\n";
  }

  ofs << std::scientific << std::setprecision(3);
  for (double E_electron = 1e12 * SI::eV; E_electron < 1e23 * SI::eV; E_electron *= 1.02) {
    ofs << E_electron / SI::eV << "   ";
    ofs << 2.0 * pow2(E_electron / SI::eV) * Q(E_electron, REDSHIFT) << "   ";
    ofs << 2.0 * pow2(E_electron / SI::eV) * Q_interpolated(E_electron, REDSHIFT) << "   ";
    ofs << 2.0 * pow2(E_electron / SI::eV) * Q_interpolated_x(E_electron, REDSHIFT) << "\n";
  }

  ofs.close();

  auto I_Q = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_electron = std::exp(lnE);
        return 2.0 * pow2(E_electron / SI::eV) * Q(E_electron, REDSHIFT);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);
  std::cout << std::setprecision(3) << "I_Q = " << I_Q << "\n";

  auto I_Q_tilde = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_electron = std::exp(lnE);
        return 2.0 * pow2(E_electron / SI::eV) * Q_interpolated(E_electron, REDSHIFT);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);
  std::cout << std::setprecision(3) << "I_Q_tilde = " << I_Q_tilde << "\n";

  auto I_Q_tilde_x = simprop::utils::QAGIntegration<double>(
      [&](double lnE) {
        auto E_electron = std::exp(lnE);
        return 2.0 * pow2(E_electron / SI::eV) * Q_interpolated_x(E_electron, REDSHIFT);
      },
      std::log(1e12 * SI::eV), std::log(1e23 * SI::eV), LIMIT, precision);
  std::cout << std::setprecision(3) << "I_Q_tilde = " << I_Q_tilde_x << "\n";

  std::cout << "written to " << filename << "\n";
}

int main() {
  try {
    simprop::utils::startup_information();

    // print_absorption_term();
    print_source_terms();

  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what();
  }
  return 0;
}