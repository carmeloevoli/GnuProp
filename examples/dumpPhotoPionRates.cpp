#include "rates.h"
#include "simprop.h"

void testInteractionRate(const std::string& filename) {
  std::string inputfile = "data/" + filename + ".bin";
  std::string outputfile = "output/" + filename + ".txt";

  gnuprop::ProductionRate gnu(inputfile);
  std::ofstream out(outputfile);
  out << "# Enu [eV] - rates[Gyr] \n";
  out << std::scientific;
  const auto units = 1. / SI::Gyr;
  auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 1000);
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    out << gnu.get(x * 1e19 * SI::eV, 1e19 * SI::eV, 0.) / units << " ";
    out << gnu.get(x * 1e20 * SI::eV, 1e20 * SI::eV, 0.) / units << " ";
    out << gnu.get(x * 1e21 * SI::eV, 1e21 * SI::eV, 0.) / units << " ";
    out << gnu.get(x * 1e22 * SI::eV, 1e22 * SI::eV, 0.) / units << " ";
    out << "\n";
  }
  LOGI << "Photo-pion interaction rates saved in " << outputfile;
}

void testInteractionRate2d(const std::string& filename) {
  std::string inputfile = "data/" + filename + ".bin";
  std::string outputfile = "output/" + filename + "_2d.txt";

  gnuprop::ProductionRate gnu(inputfile);
  std::ofstream out(outputfile);
  out << "# Enu [eV] - rates[Gyr] \n";
  out << std::scientific << std::setprecision(6);

  const double units = 1. / SI::Gyr;
  const double z = 0.;
  auto eAxis = simprop::utils::LogAxis<double>(1e17 * SI::eV, 1e23 * SI::eV, 300);
  auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 300);
  for (auto E : eAxis) {
    for (auto x : xAxis) {
      out << x << " " << E / SI::eV << " " << gnu.get(x * E, E, z) / units << "\n";
    }
  }
  LOGI << "Photo-pion interaction rates saved in " << outputfile;
}

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");

    testInteractionRate("gnuprop_photopion_neutrinos_cmb");
    testInteractionRate("gnuprop_photopion_neutrinos_ebl");
    testInteractionRate("gnuprop_photopion_gammas_cmb");
    testInteractionRate("gnuprop_photopion_pairs_cmb");
    testInteractionRate("gnuprop_photopion_gammas_ebl");
    testInteractionRate("gnuprop_photopion_pairs_ebl");

    testInteractionRate2d("gnuprop_photopion_neutrinos_cmb");
    testInteractionRate2d("gnuprop_photopion_neutrinos_ebl");

  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}