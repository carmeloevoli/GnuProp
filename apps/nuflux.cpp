#include "cosmoneutrinos.h"
#include "simprop.h"
#include "uhecr.h"

void printNuEmissivity(const beniamino::Uhecr &b, std::string filename) {
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
    beniamino::Uhecr b(std::make_unique<simprop::cosmo::Cosmology>());
    printNuEmissivity(b, "Beniamino_nuemissivity_noCutoff_2.6_0_3.txt");
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}