#include "interactions/BreitWheeler.h"
#include "rates.h"
#include "simprop.h"

void breitwheeler() {
  const auto chiAxis = simprop::utils::LogAxis<double>(1, 1e4, 1000);
  PhotonPairProduction::BreitWheeler bw;
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

void absorptionLength() {
  gnuprop::GammaAbsorptionRate gnu;
  std::ofstream out("output/gnuprop_gamma_absorption_length.txt");
  out << "# energy [TeV] - length [Mpc]\n";
  out << std::scientific;
  const auto units = SI::Mpc;
  auto eAxis = simprop::utils::LogAxis<double>(SI::TeV, 1e14 * SI::TeV, 1000);
  for (auto E : eAxis) {
    out << E / SI::TeV << "\t";
    out << (SI::cLight / gnu.get(E, 0.)) / units << " ";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();

    simprop::utils::Timer timer("main timer");

    breitwheeler();
    absorptionLength();
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}