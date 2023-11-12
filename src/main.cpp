#include "beniamino.h"
#include "cosmoneutrinos.h"
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
    std::ofstream out("output/Beniamino_characteristics_vs_energy.txt");
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

void testPhotonField() {
  simprop::photonfields::CMB cmb;
  auto epsAxis = simprop::utils::LogAxis(1e-6 * SI::eV, SI::eV, 1000);
  std::ofstream out("output/Beniamino_cmb.txt");
  const double units = 1. / SI::eV / SI::cm3;
  out << "# eps [eV] - density\n";
  for (const auto &eps : epsAxis) {
    out << std::scientific << eps / SI::eV << "\t";
    out << cmb.density(eps, 1.) / units << "\t";
    out << "\n";
  }
}

void testCrossSections() {
  KelnerAharonian2008::NeutrinoProductionSpectrum nuSpec;
  std::ofstream out("output/Beniamino_neutrino_xsecs.txt");
  out << "# x - spectrum\n";
  out << std::scientific;
  const auto units = SI::cm3 / SI::sec;
  const auto eta0 = 0.313;
  auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << nuSpec.Phi(2., x) / units << "\t";
    out << nuSpec.Phi(20., x) / units << "\t";
    out << nuSpec.nu_mu(2., x) / units << "\t";
    out << nuSpec.nu_mu(20., x) / units << "\t";
    out << nuSpec.barnu_mu(2., x) / units << "\t";
    out << nuSpec.barnu_mu(20., x) / units << "\t";
    out << nuSpec.nu_e(2., x) / units << "\t";
    out << nuSpec.nu_e(20., x) / units << "\t";
    out << nuSpec.barnu_e(2., x) / units << "\t";
    out << nuSpec.barnu_e(20., x) / units << "\t";
    out << "\n";
  }
}

void testLosses() {
  auto losses = std::make_unique<beniamino::LossesTable>();
  assert(losses->loadTable("data/SimProp_proton_losses.txt"));
  std::ofstream out("output/Beniamino_losses.txt");
  out << "# E - losses\n";
  out << std::scientific;
  const auto units = 1. / SI::year;
  auto energyAxis = simprop::utils::LogAxis<double>(1e15 * SI::eV, 1e24 * SI::eV, 1000);
  for (auto E : energyAxis) {
    out << E / SI::eV << " ";
    out << losses->beta(E) / units << " ";
    out << losses->dbdE(E) / units << " ";
    out << "\n";
  }
}

void testSpectrum(const beniamino::Beniamino &b, std::string filename) {
  auto E = simprop::utils::LogAxis(1e16 * SI::eV, 1e21 * SI::eV, 100);
  filename = "output/" + filename;
  std::ofstream out(filename);
  const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
  out << "# E [eV] - z=0 - z=1 - z=2 - z=3 - z=4 - z=5\n";
  for (const auto &E_i : E) {
    std::cout << E_i / SI::eV << "\n";
    out << std::scientific << E_i / SI::eV << "\t";
    out << b.computeFlux(E_i, 0., 1e-4) / units << "\t";
    out << b.computeFlux(E_i, 1., 1e-4) / units << "\t";
    out << b.computeFlux(E_i, 2., 1e-4) / units << "\t";
    out << b.computeFlux(E_i, 3., 1e-4) / units << "\t";
    out << b.computeFlux(E_i, 4., 1e-4) / units << "\t";
    out << b.computeFlux(E_i, 5., 1e-4) / units << "\t";
    out << "\n";
  }
}

void printNuEmissivity(const beniamino::Beniamino &b, std::string filename) {
  auto E = simprop::utils::LogAxis(1e17 * SI::eV, 1e21 * SI::eV, 4 * 16);
  beniamino::CosmoNeutrinos cosmonus(b);
  filename = "output/" + filename;
  std::ofstream out(filename);
  const double units = 1. / SI::eV / SI::m3 / SI::sec;
  out << "# E [eV] - z=0\n";
  for (const auto &E_i : E) {
    std::cout << std::scientific << E_i / SI::eV << "\n";
    out << std::scientific << E_i / SI::eV << "\t";
    out << cosmonus.nuEmissivity(E_i, 0.1) / units << "\t";
    out << cosmonus.nuEmissivity(E_i, 0.5) / units << "\t";
    out << cosmonus.nuEmissivity(E_i, 1.0) / units << "\t";
    out << cosmonus.nuEmissivity(E_i, 2.0) / units << "\t";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    beniamino::Beniamino b(std::make_unique<simprop::cosmo::Cosmology>());
    // testCharacteristics(b);
    // testJacobian(b);
    // testPhotonField();
    testCrossSections();
    // testSpectrum(b, "Beniamino_spectrum_noCutoff_2.6_0_3.txt");
    // printNuEmissivity(b, "Beniamino_nuemissivity_noCutoff_2.6_0_3.txt");
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}