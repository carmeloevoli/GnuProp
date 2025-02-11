#include "rates.h"
#include "simprop.h"

void absorptionRate(const std::string& inputfile, const std::string& outputfile) {
  gnuprop::AbsorptionRate rate("data/" + inputfile);
  std::ofstream out("output/" + outputfile);
  out << "# energy [eV] - absorption rate [Gyr-1]\n";
  out << std::scientific;
  const auto units = 1. / SI::Gyr;
  auto eAxis = simprop::utils::LogAxis<double>(1e10 * SI::eV, 1e20 * SI::eV, 100);
  for (auto E : eAxis) {
    out << E / SI::eV << "\t";
    out << rate.get(E, 0.) / units << " ";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");

    // rates
    absorptionRate("gnuprop_absorption_gammas_cmb.bin", "gnuprop_absorption_gammas_cmb.txt");
    absorptionRate("gnuprop_absorption_gammas_ebl.bin", "gnuprop_absorption_gammas_ebl.txt");
    absorptionRate("gnuprop_absorption_pairs_cmb.bin", "gnuprop_absorption_pairs_cmb.txt");
    absorptionRate("gnuprop_absorption_pairs_ebl.bin", "gnuprop_absorption_pairs_ebl.txt");

  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}