#include "interactions/InverseCompton.h"
#include "interactions/PhotoPair.h"
#include "simprop.h"

// #include "interactions/InverseCompton.h"
// #include "interactions/PhotoPair.h"

// void testCrossSections() {
//   {
//     interactions::PhotoPionNeutrinos nuSpec;
//     std::string filename = "output/gnuprop_neutrino_KA_phi.txt";
//     std::ofstream out(filename);
//     out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
//     out << std::scientific;
//     const auto units = SI::cm3 / SI::sec;
//     auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
//     for (auto x : xAxis) {
//       out << std::scientific << x << "\t";
//       out << nuSpec.Phi(0.5, x) / units << "\t";
//       out << nuSpec.Phi(2., x) / units << "\t";
//       out << nuSpec.Phi(20., x) / units << "\t";
//       out << nuSpec.Phi(30., x) / units << "\t";
//       out << "\n";
//     }
//     LOGI << "Neutrino cross-sections saved in " << filename;
//   }
//   {
//     interactions::PhotoPionGammas gammaSpec;
//     std::string filename = "output/gnuprop_gamma_KA_phi.txt";
//     std::ofstream out(filename);
//     out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
//     out << std::scientific;
//     const auto units = SI::cm3 / SI::sec;
//     auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
//     for (auto x : xAxis) {
//       out << std::scientific << x << "\t";
//       out << gammaSpec.Phi(0.5, x) / units << "\t";
//       out << gammaSpec.Phi(2., x) / units << "\t";
//       out << gammaSpec.Phi(20., x) / units << "\t";
//       out << gammaSpec.Phi(30., x) / units << "\t";
//       out << "\n";
//     }
//     LOGI << "Gamma cross-sections saved in " << filename;
//   }
//   {
//     interactions::PhotoPionPairs pairSpec;
//     std::string filename = "output/gnuprop_pair_KA_phi.txt";
//     std::ofstream out(filename);
//     out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
//     out << std::scientific;
//     const auto units = SI::cm3 / SI::sec;
//     auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
//     for (auto x : xAxis) {
//       out << std::scientific << x << "\t";
//       out << pairSpec.Phi(0.5, x) / units << "\t";
//       out << pairSpec.Phi(2., x) / units << "\t";
//       out << pairSpec.Phi(20., x) / units << "\t";
//       out << pairSpec.Phi(30., x) / units << "\t";
//       out << "\n";
//     }
//     LOGI << "Pair cross-sections saved in " << filename;
//   }
// }

void inversecompton() {
  const auto chiAxis = simprop::utils::LogAxis<double>(0.1, 1e4, 1000);
  Interactions::InverseCompton ic;
  std::ofstream out("output/gnuprop_xsecs_inversecompton.txt");
  out << "# chi - s [GeV^2] - sigma [mbarn]\n";
  out << std::scientific;
  for (auto chi : chiAxis) {
    auto s = chi * pow2(SI::electronMassC2);
    out << chi << "\t";
    out << s / SI::GeV2 << "\t";
    out << ic.sigma_com(s) / SI::mbarn << "\t";
    out << "\n";
  }
}

void breitwheeler() {
  const auto chiAxis = simprop::utils::LogAxis<double>(0.1, 1e4, 1000);
  Interactions::PhotoPair photoPair;
  std::ofstream out("output/gnuprop_xsecs_breitwheeler.txt");
  out << "# chi - s [GeV^2] - sigma [mbarn]\n";
  out << std::scientific;
  for (auto chi : chiAxis) {
    auto s = chi * pow2(2. * SI::electronMassC2);
    out << chi << "\t";
    out << s / SI::GeV2 << "\t";
    out << photoPair.sigma_com(s) / SI::mbarn << "\t";
    out << "\n";
  }
}

void inversecompton_differential() {
  const auto eElectron = 1e16 * SI::eV;
  const auto eBkg = 2.35e-4 * SI::eV;
  const auto eGammaAxis = simprop::utils::LogAxis<double>(1e12 * SI::eV, 1e16 * SI::eV, 1000);
  const auto units = SI::mbarn;
  Interactions::InverseCompton ic;
  std::ofstream out("output/gnuprop_xsecs_inversecompton_differential.txt");
  out << "# E_gamma [eV] - E dsigmadE [mbarn]\n";
  out << std::scientific;
  for (auto eGamma : eGammaAxis) {
    out << eGamma / SI::eV << "\t";
    out << eGamma * ic.dsigma_dE(eElectron, eBkg, eGamma) / units << "\t";
    out << "\n";
  }
}

void breitwheeler_differential() {
  const auto eElectron = 1e4 * SI::TeV;
  const auto eBkg = 2.35e-4 * SI::eV;
  const auto eGammaAxis = simprop::utils::LogAxis<double>(1e4 * SI::TeV, 1e9 * SI::TeV, 1000);
  const auto units = SI::mbarn;
  Interactions::PhotoPair photoPair;
  std::ofstream out("output/gnuprop_xsecs_breitwheeler_differential.txt");
  out << "# E_gamma [eV] - E dsigmadE [mbarn]\n";
  out << std::scientific;
  for (auto eGamma : eGammaAxis) {
    out << eGamma / SI::eV << "\t";
    out << eElectron * photoPair.dsigma_dE(eGamma, eBkg, eElectron) / units << "\t";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");

    // cross-sections
    breitwheeler();
    breitwheeler_differential();
    inversecompton();
    inversecompton_differential();
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}