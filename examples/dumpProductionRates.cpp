#include "rates.h"
#include "simprop.h"

void dumpRates(const std::string& filename) {
  std::string inputfile = "tables/" + filename + ".bin";
  std::string outputfile = "output/" + filename + ".txt";

  gnuprop::ProductionRate gnu(inputfile);
  std::ofstream out(outputfile);
  out << "# x - rates [Gyr^-1]\n";
  out << std::scientific << std::setprecision(6);
  const auto units = 1. / SI::Gyr;
  auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 1000);
  std::vector<double> eAxis{1e10, 1e12, 1e14, 1e16, 1e18, 1e20, 1e22, 1e24};
  for (auto x : xAxis) {
    out << std::scientific << x << "\t";
    for (auto E : eAxis) out << gnu.get(x * E * SI::eV, E * SI::eV, 0.) / units << "\t";
    out << "\n";
  }
  LOGI << "Interaction rates saved in " << outputfile;
}

void dumpRates2D(const std::string& filename) {
  std::string inputfile = "tables/" + filename + ".bin";
  std::string outputfile = "output/" + filename + "_2d.txt";

  gnuprop::ProductionRate gnu(inputfile);
  std::ofstream out(outputfile);
  out << "# rates [Gyr^-1]\n";
  out << std::scientific << std::setprecision(6);

  const double units = 1. / SI::Gyr;
  const double z = 0.;
  auto eAxis = simprop::utils::LogAxis<double>(1e12 * SI::eV, 1e24 * SI::eV, 800);
  auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 400);
  for (auto E : eAxis) {
    for (auto x : xAxis) {
      out << x << " " << E / SI::eV << " " << gnu.get(x * E, E, z) / units << "\n";
    }
  }
  LOGI << "Interaction rates saved in " << outputfile;
}
int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");

    dumpRates("gnuprop_photopair_pairs_cmb");
    dumpRates("gnuprop_photopair_pairs_ebl");
    dumpRates("gnuprop_photopion_neutrinos_cmb");
    dumpRates("gnuprop_photopion_neutrinos_ebl");
    dumpRates("gnuprop_photopion_gammas_cmb");
    dumpRates("gnuprop_photopion_gammas_ebl");
    dumpRates("gnuprop_photopion_pairs_cmb");
    dumpRates("gnuprop_photopion_pairs_ebl");
    dumpRates("gnuprop_inversecompton_gammas_cmb");
    dumpRates("gnuprop_inversecompton_gammas_ebl");
    dumpRates("gnuprop_inversecompton_pairs_cmb");
    dumpRates("gnuprop_inversecompton_pairs_ebl");

    dumpRates2D("gnuprop_photopair_pairs_cmb");
    dumpRates2D("gnuprop_photopair_pairs_ebl");
    dumpRates2D("gnuprop_inversecompton_gammas_cmb");
    dumpRates2D("gnuprop_inversecompton_gammas_ebl");
    dumpRates2D("gnuprop_inversecompton_pairs_cmb");
    dumpRates2D("gnuprop_inversecompton_pairs_ebl");

  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}