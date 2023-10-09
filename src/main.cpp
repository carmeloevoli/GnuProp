#include <fstream>

#include "beniamino.h"
#include "simprop/utils/numeric.h"

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
int main() {
  try {
    printSpectrum(beniamino::Beniamino({2.6, 0., -1., 3.}), "Beniamino_spectrum_2.6_0.txt");
    printSpectrum(beniamino::Beniamino({2.6, 3., -1., 3.}), "Beniamino_spectrum_2.6_3.txt");
    printSpectrum(beniamino::Beniamino({2.6, -3., -1., 3.}), "Beniamino_spectrum_2.6_-3.txt");
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}