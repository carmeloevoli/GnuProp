#include "KelnerAharonian2008.h"
#include "simprop.h"

void testCrossSections() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  std::ofstream out("output/Beniamino_neutrino_xsecs.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << nuSpec.Phi(0.5, x) / units << "\t";
    out << nuSpec.Phi(2., x) / units << "\t";
    out << nuSpec.Phi(20., x) / units << "\t";
    out << nuSpec.Phi(40., x) / units << "\t";
    out << "\n";
  }
}

void testCrossSectionMaps() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  std::ofstream out("output/Beniamino_neutrino_xsecs_map.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  const auto epsCmb = 6.3e-4 * SI::eV;
  const auto epsIr = 1e-2 * SI::eV;
  auto EpAxis = simprop::utils::LogAxis<double>(1e16 * SI::eV, 1e22 * SI::eV, 251);
  auto EnuAxis = simprop::utils::LogAxis<double>(1e15 * SI::eV, 1e22 * SI::eV, 251);
  for (auto Enu : EnuAxis) {
    for (auto Ep : EpAxis) {
      auto x = Enu / Ep;
      {
        auto eta = 4. * epsCmb * Ep / pow2(SI::protonMassC2);
        out << std::scientific << x << " " << eta << " ";
        out << x * x * nuSpec.Phi(eta, x) / units << " ";
      }
      {
        auto eta = 4. * epsIr * Ep / pow2(SI::protonMassC2);
        out << x * x * nuSpec.Phi(eta, x) / units << " ";
      }
      out << "\n";
    }
  }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    testCrossSections();
    testCrossSectionMaps();
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}
