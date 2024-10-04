#include "losses.h"
#include "simprop.h"

void testLosses() {
  beniamino::EnergyLosses losses;
  std::ofstream out("output/GammaProp_proton_losses.txt");
  auto energyAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e24 * SI::eV, 1000);
  const auto units = SI::Mpc;
  for (const auto &E : energyAxis) {
    out << std::scientific << E / SI::eV << "\t";
    out << SI::cLight / losses.beta(E) / units << "\t";
    out << losses.beta(E) * SI::Gyr << "\t";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    testLosses();
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }

  return EXIT_SUCCESS;
}