#include "beniamino.h"

int main() {
  try {
    beniamino::Beniamino b({2.7, 0., -1., 1.});

    auto E = utils::LogAxis(1e17 * SI::eV, 1e21 * SI::eV, 16 * 4);
    std::ofstream out("output/SimProp_spectrum.txt");
    const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
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
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}