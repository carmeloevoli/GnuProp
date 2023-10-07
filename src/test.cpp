#include <fstream>
#include <iostream>

#include "beniamino.h"

void testCharacteristics(const beniamino::Beniamino &b) {
  {
    std::ofstream out("output/Beniamino_characteristics_vs_redshift.txt");
    auto z = utils::LogAxis(1e-4, 5., 1000);
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
    auto E = utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
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
    auto z = utils::LogAxis(1e-4, 1e1, 1000);
    out << "# z\n";
    for (const auto &z_i : z) {
      std::cout << z_i << "\n";
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
    auto E = utils::LogAxis(1e16 * SI::eV, 1e22 * SI::eV, 1000);
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

int main() {
  std::cout << "Hello Beniamino!\n";
  beniamino::Beniamino b;
  testCharacteristics(b);
  testJacobian(b);
}