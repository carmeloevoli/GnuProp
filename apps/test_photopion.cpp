#include <fstream>

#include "gnuprop.h"
#include "interactions/PhotoPion.h"
#include "simprop.h"

#define LIMIT 10000
#define REDSHIFT 4.5

const auto m_phField = std::make_unique<simprop::photonfields::CMB>();
const auto precision = 1e-3;
const auto m_photoPion = std::make_unique<interactions::PhotoPionNeutrinos>();
const auto q_nu =
    std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_neutrinos_cmb.bin");

#undef GAUSSIAN

#ifdef GAUSSIAN
double f_p(double E_gamma) {
  const double E_0 = 1e16 * SI::eV;
  const double sigma = 0.2 * E_0;
  return 1e20 * std::exp(-pow2((E_gamma - E_0) / sigma));
}
#else
double f_p(double E_p) {
  const double N_0 = 1e20;
  const double EeV = 1e18 * SI::eV;
  return N_0 * std::pow(E_p / SI::PeV, -1.0) * std::exp(-E_p / EeV);
}
#endif

double R(double E_nu, double E_p, double z = 0.) {
  if (E_nu > E_p) return 0.;

  const auto m_pi = SI::pionMassC2;
  const auto m_p = SI::protonMassC2;
  const auto epsTh = (pow2(m_pi) + 2. * m_p * m_pi) / 4. / E_p;
  const auto epsMin = std::max(epsTh, m_phField->getMinPhotonEnergy());
  const auto epsMax = m_phField->getMaxPhotonEnergy();

  if (epsMin > epsMax) return 0.;

  auto integrand = [&](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * E_p / pow2(SI::protonMassC2);
    auto value = m_phField->density(epsilon, z) * m_photoPion->Phi(eta, E_nu / E_p);
    return epsilon * value;
  };

  const auto a = std::log(epsMin);
  const auto b = std::log(epsMax);
  const size_t N = 10000;
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, precision);

  return std::max(value, 0.);
}

double Q(double E_nu, double z = 0.) {
  auto integrand = [&](double lnEp) {
    const auto E_p = std::exp(lnEp);
    return R(E_nu, E_p, z) * f_p(E_p);
  };

  const auto a = std::log(E_nu);
  const auto b = std::log(1e22 * SI::eV);
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, LIMIT, precision);

  return std::max(value, 0.);
}

double Q_interpolated(double E_nu, double z = 0.) {
  auto integrand = [&](double lnEp) {
    const auto E_p = std::exp(lnEp);
    return q_nu->get(E_nu / E_p, E_p, z) * f_p(E_p);
  };

  const auto a = std::log(E_nu);
  const auto b = std::log(1e22 * SI::eV);
  auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, LIMIT, 1e-2);

  return std::max(value, 0.);
}

void print_source_terms(std::string filename = "output/test_photopion_source.txt") {
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cerr << "Error: could not open output.txt for writing\n";
  }

  ofs << std::scientific << std::setprecision(3);
  for (double E_nu = 1e12 * SI::eV; E_nu < 1e23 * SI::eV; E_nu *= 1.1) {
    ofs << E_nu / SI::eV << "   ";
    LOGI << E_nu / SI::eV;
    ofs << pow2(E_nu / SI::eV) * Q(E_nu, REDSHIFT) << "   ";
    ofs << pow2(E_nu / SI::eV) * Q_interpolated(E_nu, REDSHIFT) << "\n";
  }

  ofs.close();

  std::cout << "written to " << filename << "\n";
}

int main() {
  try {
    simprop::utils::startup_information();

    print_source_terms();

  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what();
  }
  return 0;
}