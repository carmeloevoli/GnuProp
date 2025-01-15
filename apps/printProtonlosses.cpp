#include "losses.h"
#include "simprop.h"

void testLosses() {
  gnuprop::EnergyLosses losses;

  // Open output file for writing
  const std::string outputFile = "output/gnuprop_proton_losses.txt";
  std::ofstream out(outputFile);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output file: " + outputFile);
  }

  // Configure output format
  out << "# Energy (eV)\tMean Free Path (Mpc)\tBeta (Gyr^-1)\n";
  out << std::scientific << std::setprecision(6);

  // Generate energy axis
  const auto energyAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e24 * SI::eV, 1000);
  const double units = SI::Mpc;

  // Compute and write results
  for (const auto& energy : energyAxis) {
    const double meanFreePath = SI::cLight / losses.beta(energy) / units;
    const double betaValue = losses.beta(energy) * SI::Gyr;

    out << energy / SI::eV << "\t" << meanFreePath << "\t" << betaValue << "\n";
  }

  LOGI << "Results written to " << outputFile;
}

int main() {
  try {
    // Display startup information
    simprop::utils::startup_information();

    // Measure execution time
    simprop::utils::Timer timer("main timer");

    // Run the test
    testLosses();
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
