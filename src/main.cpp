#include "beniamino.h"
#include "simprop.h"

void testCharacteristics(const beniamino::Beniamino &b) {
  {
    std::ofstream out("output/Beniamino_characteristics_vs_redshift.txt");
    auto z = simprop::utils::LogAxis(1e-4, 5., 1000);
    out << "# z - 10^17 eV - 10^18 eV - 10^19 eV - 10^20 eV - 10^21 eV\n";
    for (const auto &z_i : z) {
      out << std::scientific << z_i << "\t";
      out << b.generationEnergy(1e17 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e18 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e19 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e20 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(1e21 * SI::eV, 0., z_i, 1e-5) / SI::eV << "\t";
      out << "\n";
    }
    out.close();
  }
  {
    std::ofstream out("output/Beniamino_characteristics_vs_energy.txt");
    auto E = simprop::utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
    out << "# E - 0.05 - 0.5 - 1 - 2 - 3 - 5\n";
    for (const auto &E_i : E) {
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 0.05, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 0.5, 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 1., 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 2., 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 3., 1e-5) / SI::eV << "\t";
      out << b.generationEnergy(E_i, 0., 5., 1e-5) / SI::eV << "\t";
      out << "\n";
    }
    out.close();
  }
}

void testJacobian(const beniamino::Beniamino &b) {
  {
    std::ofstream out("output/Beniamino_jacobian_vs_redshift.txt");
    auto z = simprop::utils::LogAxis(1e-4, 1e1, 1000);
    out << "# z - 10^17 eV - 10^18 eV - 10^19 eV - 10^20 eV - 10^21 eV\n";
    for (const auto &z_i : z) {
      out << std::scientific << z_i << "\t";
      out << b.dilationFactor(1e17 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e18 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e19 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e20 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << b.dilationFactor(1e21 * SI::eV, 0., z_i, 1e-6) << "\t";
      out << "\n";
    }
    out.close();
  }
  {
    std::ofstream out("output/Beniamino_jacobian_vs_energy.txt");
    auto E = simprop::utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
    out << "# E - 0.05 - 0.5 - 1 - 2 - 3 - 5\n";
    for (const auto &E_i : E) {
      std::cout << E_i / SI::eV << "\n";
      out << std::scientific << E_i / SI::eV << "\t";
      out << b.dilationFactor(E_i, 0., 0.05, 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 0.5, 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 1., 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 2., 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 3., 1e-6) << "\t";
      out << b.dilationFactor(E_i, 0., 5., 1e-6) << "\t";
      out << "\n";
    }
    out.close();
  }
}

void printSpectrum(const beniamino::Beniamino &b, std::string filename) {
  auto E = simprop::utils::LogAxis(1e17 * SI::eV, 1e21 * SI::eV, 16 * 4);
  filename = "output/" + filename;
  std::ofstream out(filename);
  const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
  out << "# E [eV] - z=0 - z=1 - z=2 - z=3 - z=4 - z=5\n";
  for (const auto &E_i : E) {
    std::cout << E_i / SI::eV << "\n";
    out << std::scientific << E_i / SI::eV << "\t";
    out << b.computeFlux(E_i, 0.) / units << "\t";
    out << b.computeFlux(E_i, 1.) / units << "\t";
    out << b.computeFlux(E_i, 2.) / units << "\t";
    out << b.computeFlux(E_i, 3.) / units << "\t";
    out << b.computeFlux(E_i, 4.) / units << "\t";
    out << b.computeFlux(E_i, 5.) / units << "\t";
    out << "\n";
  }
}

void printNuSpectrum(const beniamino::Beniamino &b, std::string filename) {
  auto E = simprop::utils::LogAxis(1e17 * SI::eV, 1e18 * SI::eV, 3);
  // beniamino::CosmoNeutrinos nus(b);
  // filename = "output/" + filename;
  // std::ofstream out(filename);
  // const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
  // out << "# E [eV] - z=0\n";
  // for (const auto &E_i : E) {
  //   std::cout << E_i / SI::eV << "\n";
  //   out << std::scientific << E_i / SI::eV << "\t";
  //   out << nus.computeNeutrinoFlux(E_i, 1.) / units << "\t";
  //   out << "\n";
  // }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    beniamino::Beniamino b(std::make_unique<simprop::cosmo::Cosmology>());
    testCharacteristics(b);
    testJacobian(b);
    printSpectrum(b, "Beniamino_spectrum_noCutoff_2.6_0_3.txt");

    // printSpectrum(beniamino::Beniamino({2.6, 0., -1., 3.}), "Beniamino_spectrum_2.6_0.txt");
    // printSpectrum(beniamino::Beniamino({2.6, 3., -1., 3.}), "Beniamino_spectrum_2.6_3.txt");
    // printSpectrum(beniamino::Beniamino({2.6, -3., -1., 3.}), "Beniamino_spectrum_2.6_-3.txt");
    // printNuSpectrum(beniamino::Beniamino({2.6, 0., -1., 3.}), "Beniamino_nuspectrum_2.6_0.txt");

  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}