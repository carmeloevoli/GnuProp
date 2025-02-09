#include "interactions/GammaPairProduction.h"
#include "rates.h"
#include "simprop.h"

void breitwheeler() {
  const auto chiAxis = simprop::utils::LogAxis<double>(1, 1e4, 1000);
  Interactions::GammaPairProduction bw;
  std::ofstream out("output/gnuprop_xsecs_breitwheeler.txt");
  out << "# chi - s [GeV^2] - sigma [mbarn]\n";
  out << std::scientific;
  for (auto chi : chiAxis) {
    auto s = chi * pow2(2. * SI::electronMassC2);
    out << chi << "\t";
    out << s / SI::GeV2 << "\t";
    out << bw.sigmaInCoMFrame(s) / SI::mbarn << "\t";
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

    breitwheeler();
    absorptionRate("data/gnuprop_absorption_gammas_cmb.bin",
                   "output/gnuprop_gamma_absorption_cmb.txt");
    absorptionRate("data/gnuprop_absorption_gammas_ebl.bin",
                   "output/gnuprop_gamma_absorption_ebl.txt");

  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}