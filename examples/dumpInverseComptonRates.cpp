#include "rates.h"
#include "simprop.h"

void testInteractionRate(const std::string& filename) {
  std::string inputfile = "data/" + filename + ".bin";
  std::string outputfile = "output/" + filename + ".txt";

  gnuprop::ProductionRate gnu(inputfile);
  std::ofstream out(outputfile);
  out << "# E [eV] - rates[Gyr] \n";
  out << std::scientific;
  const auto units = 1. / SI::Gyr;
  auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << gnu.get(x * 1e12 * SI::eV, 1e15 * SI::eV, 0.) / units << " ";
    out << gnu.get(x * 1e13 * SI::eV, 1e16 * SI::eV, 0.) / units << " ";
    out << gnu.get(x * 1e14 * SI::eV, 1e17 * SI::eV, 0.) / units << " ";
    out << gnu.get(x * 1e15 * SI::eV, 1e18 * SI::eV, 0.) / units << " ";
    out << "\n";
  }
  LOGI << "Photo-pair interaction rates saved in " << outputfile;
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");

    testInteractionRate("gnuprop_IC_gammas_cmb");
    testInteractionRate("gnuprop_IC_gammas_ebl");
    testInteractionRate("gnuprop_IC_pairs_cmb");
    testInteractionRate("gnuprop_IC_pairs_ebl");

  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}