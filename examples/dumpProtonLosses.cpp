#include "rates.h"
#include "simprop.h"

using cmb = simprop::photonfields::CMB;

void protonLosses() {
  simprop::losses::PairProductionLosses losses_pair(std::make_unique<cmb>());
  simprop::losses::PhotoPionContinuousLosses losses_photopi(std::make_unique<cmb>());

  // Open output file for writing
  const std::string outputFile = "output/gnuprop_proton_losses_sirente.txt";
  std::ofstream out(outputFile);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output file: " + outputFile);
  }

  // Configure output format
  out << "# Energy (eV)\tpair-production (Gyr^-1)\tphoto-pion (Gyr^-1)\n";
  out << std::scientific << std::setprecision(6);

  // Generate energy axis
  const auto energyAxis = simprop::utils::LogAxis(1e17 * SI::eV, 1e24 * SI::eV, 300);
  const auto units = 1. / SI::Gyr;

  // Compute and write results
  for (const auto& energy : energyAxis) {
    const auto gammap = energy / SI::protonMassC2;
    auto beta_ee = losses_pair.beta(simprop::proton, gammap) / units;
    auto beta_pi = losses_photopi.beta(simprop::proton, gammap) / units;
    out << energy / SI::eV << "\t" << beta_ee << "\t" << beta_pi << "\n";
  }

  LOGI << "Results written to " << outputFile;
}

void protonLossesInterpolated() {
  gnuprop::ProtonLossRate losses_pair("data/gnuprop_proton_losses_pair.bin");
  gnuprop::ProtonLossRate losses_pi("data/gnuprop_proton_losses_photopion.bin");

  // Open output file for writing
  const std::string outputFile = "output/gnuprop_proton_losses_interpolated.txt";
  std::ofstream out(outputFile);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output file: " + outputFile);
  }

  // Configure output format
  out << "# Energy (eV)\tbeta (Gyr^-1)\n";
  out << std::scientific << std::setprecision(6);

  // Generate energy axis
  const auto energyAxis = simprop::utils::LogAxis(1e17 * SI::eV, 1e24 * SI::eV, 10000);
  const auto units = 1. / SI::Gyr;

  // Compute and write results
  for (const auto& energy : energyAxis) {
    auto beta = (losses_pair.beta(energy) + losses_pi.beta(energy)) / units;
    out << energy / SI::eV << "\t" << beta << "\n";
  }

  LOGI << "Results written to " << outputFile;
}

int main() {
  try {
    // Display startup information
    simprop::utils::startup_information();

    // Timer
    simprop::utils::Timer timer("proton losses");

    // Compute the tables
    protonLosses();
    protonLossesInterpolated();
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
