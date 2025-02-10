#include "interactions/InverseCompton.h"
#include "rates.h"
#include "simprop.h"

void inversecompton() {
  const auto chiAxis = simprop::utils::LogAxis<double>(0.1, 1e4, 1000);
  Interactions::InverseCompton ic;
  std::ofstream out("output/gnuprop_xsecs_inversecompton.txt");
  out << "# chi - s [GeV^2] - sigma [mbarn]\n";
  out << std::scientific;
  for (auto chi : chiAxis) {
    auto s = chi * pow2(SI::electronMassC2);
    out << chi << "\t";
    out << s / SI::GeV2 << "\t";
    out << ic.sigmaInCoMFrame(s) / SI::mbarn << "\t";
    out << "\n";
  }
}

void absorptionRate(const std::string& inputfile, const std::string& outputfile) {
  gnuprop::GammaAbsorptionRate rate(inputfile);
  std::ofstream out(outputfile);
  out << "# energy [eV] - absorption rate [Gyr-1]\n";
  out << std::scientific;
  const auto units = 1. / SI::Gyr;
  auto eAxis = simprop::utils::LogAxis<double>(1e12 * SI::eV, 1e24 * SI::eV, 1000);
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

    inversecompton();
    absorptionRate("data/gnuprop_absorption_pairs_cmb.bin",
                   "output/gnuprop_pair_absorption_cmb.txt");
    absorptionRate("data/gnuprop_absorption_pairs_ebl.bin",
                   "output/gnuprop_pair_absorption_ebl.txt");
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}