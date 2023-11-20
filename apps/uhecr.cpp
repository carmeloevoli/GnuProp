#include "uhecr.h"

#include "simprop.h"

void testSpectrum(const beniamino::Uhecr &b, std::string filename) {
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

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    beniamino::Uhecr b(std::make_unique<simprop::cosmo::Cosmology>());
    testSpectrum(b, "Beniamino_spectrum_noCutoff_2.6_0_3.txt");
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}