#include "gnuprop.h"
#include "interactions/KelnerAharonian2008.h"
#include "rates.h"
#include "simprop.h"

void testCrossSections() {
  {
    KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
    std::string filename = "output/gnuprop_neutrino_KA_phi.txt";
    std::ofstream out(filename);
    out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
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
    LOGI << "Neutrino cross-sections saved in " << filename;
  }
  {
    KelnerAharonian2008::GammaSpectrum gammaSpec;
    std::string filename = "output/gnuprop_gamma_KA_phi.txt";
    std::ofstream out(filename);
    out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
    out << std::scientific;
    const auto units = SI::cm3 / SI::sec;
    auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
    for (auto x : xAxis) {
      out << std::scientific << x << "\t";
      out << gammaSpec.Phi(0.5, x) / units << "\t";
      out << gammaSpec.Phi(2., x) / units << "\t";
      out << gammaSpec.Phi(20., x) / units << "\t";
      out << gammaSpec.Phi(30., x) / units << "\t";
      out << "\n";
    }
    LOGI << "Gamma cross-sections saved in " << filename;
  }
  {
    KelnerAharonian2008::PairSpectrum pairSpec;
    std::string filename = "output/gnuprop_pair_KA_phi.txt";
    std::ofstream out(filename);
    out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
    out << std::scientific;
    const auto units = SI::cm3 / SI::sec;
    auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
    for (auto x : xAxis) {
      out << std::scientific << x << "\t";
      out << pairSpec.Phi(0.5, x) / units << "\t";
      out << pairSpec.Phi(2., x) / units << "\t";
      out << pairSpec.Phi(20., x) / units << "\t";
      out << pairSpec.Phi(30., x) / units << "\t";
      out << "\n";
    }
    LOGI << "Pair cross-sections saved in " << filename;
  }
  {
    KelnerAharonian2008::ElectronSpectrum electronSpec;
    std::string filename = "output/gnuprop_electron_KA_phi.txt";
    std::ofstream out(filename);
    out << "# x | phi (0.5) | phi (2) | phi(20) | phi(30) [cm3/s]\n";
    out << std::scientific;
    const auto units = SI::cm3 / SI::sec;
    auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
    for (auto x : xAxis) {
      out << std::scientific << x << "\t";
      out << electronSpec.Phi(0.5, x) / units << "\t";
      out << electronSpec.Phi(2., x) / units << "\t";
      out << electronSpec.Phi(20., x) / units << "\t";
      out << electronSpec.Phi(30., x) / units << "\t";
      out << "\n";
    }
    LOGI << "Pair cross-sections saved in " << filename;
  }
}

void testPartialCrossSections() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  std::ofstream out("output/gnuprop_neutrino_xsecs_partial.txt");
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
  gnuprop::NeutrinoProductionRate gnu;
  {
    std::ofstream out("output/gnuprop_neutrino_rate_xaxis.txt");
    out << "# Enu [eV] - rates[Gyr] \n";
    out << std::scientific;
    const auto units = 1. / SI::Gyr;
    auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
    for (auto x : xAxis) {
      out << std::scientific << x << "\t";
      out << gnu.get(x * 1e19 * SI::eV, 1e19 * SI::eV, 0.) / units << " ";
      out << gnu.get(x * 1e20 * SI::eV, 1e20 * SI::eV, 0.) / units << " ";
      out << gnu.get(x * 1e21 * SI::eV, 1e21 * SI::eV, 0.) / units << " ";
      out << gnu.get(x * 1e22 * SI::eV, 1e22 * SI::eV, 0.) / units << " ";
      out << "\n";
    }
  }
  {
    std::ofstream out("output/gnuprop_neutrino_rate_zaxis.txt");
    out << "# Enu [eV] - rates[Gyr] \n";
    out << std::scientific;
    const auto units = 1. / SI::Gyr;
    const auto x = 0.1;
    auto zAxis = simprop::utils::LinAxis<double>(0, 10, 100);
    for (auto z : zAxis) {
      out << std::scientific << z << "\t";
      out << gnu.get(x * 1e19 * SI::eV, 1e19 * SI::eV, z) / units << " ";
      out << gnu.get(x * 1e20 * SI::eV, 1e20 * SI::eV, z) / units << " ";
      out << gnu.get(x * 1e21 * SI::eV, 1e21 * SI::eV, z) / units << " ";
      out << gnu.get(x * 1e22 * SI::eV, 1e22 * SI::eV, z) / units << " ";
      out << "\n";
    }
  }
  {
    std::ofstream out("output/gnuprop_neutrino_rate_eaxis.txt");
    out << "# Enu [eV] - rates[Gyr] \n";
    out << std::scientific;
    const auto units = 1. / SI::Gyr;
    const auto z = 0;
    auto eAxis = simprop::utils::LogAxis<double>(1e17 * SI::eV, 1e22 * SI::eV, 1000);
    for (auto E : eAxis) {
      out << std::scientific << E / SI::eV << "\t";
      out << gnu.get(1e-1 * E, E, z) / units << " ";
      out << gnu.get(1e-2 * E, E, z) / units << " ";
      out << gnu.get(1e-3 * E, E, z) / units << " ";
      out << gnu.get(1e-4 * E, E, z) / units << " ";
      out << "\n";
    }
  }
}

int main() {
  try {
    simprop::utils::startup_information();

    simprop::utils::Timer timer("main timer");

    testCrossSections();
    testPartialCrossSections();
    // testInteractionRate();
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}