#include "KelnerAharonian2008.h"
#include "gammaprop.h"
#include "simprop.h"

void testCrossSections() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  std::ofstream out("output/GammaProp_neutrino_xsecs.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << nuSpec.Phi(0.5, x) / units << "\t";
    out << nuSpec.Phi(2., x) / units << "\t";
    out << nuSpec.Phi(20., x) / units << "\t";
    out << nuSpec.Phi(30., x) / units << "\t";
    out << "\n";
  }
}

void testPartialCrossSections() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  std::ofstream out("output/GammaProp_neutrino_xsecs_partial.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  const auto eta = 2.;
  auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << nuSpec.nu_e(eta, x) / units << "\t";
    out << nuSpec.nu_mu(eta, x) / units << "\t";
    out << nuSpec.barnu_e(eta, x) / units << "\t";
    out << nuSpec.barnu_mu(eta, x) / units << "\t";
    out << "\n";
  }
}

void testInteractionRate() {
  gammaprop::GammaProp g(std::make_unique<simprop::cosmo::Cosmology>());
  std::ofstream out("output/GammaProp_neutrino_rate.txt");
  out << "# Enu [eV] - rates[Gyr] \n";
  out << std::scientific;
  const auto units = 1. / SI::Gyr;
  auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << g.nuInteractionRate(x * 1e19 * SI::eV, 1e19 * SI::eV, 0.) / units << " ";
    out << g.nuInteractionRate(x * 1e20 * SI::eV, 1e20 * SI::eV, 0.) / units << " ";
    out << g.nuInteractionRate(x * 1e21 * SI::eV, 1e21 * SI::eV, 0.) / units << " ";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    testCrossSections();
    testPartialCrossSections();
    testInteractionRate();
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}
