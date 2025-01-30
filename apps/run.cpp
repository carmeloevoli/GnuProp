#include "gnuprop.h"
#include "simprop.h"

int main() {
  try {
    simprop::utils::startup_information();
    {
      simprop::utils::Timer timer("timer 0");
      gnuprop::GnuProp g(std::make_unique<simprop::cosmo::Cosmology>());
      g.setRedshiftSize(1e3);
      g.setRedshiftMax(5.0);
      g.disablePhotoPion();
      g.build();
      g.evolve(0.);
      g.dump("gnuprop_spectrum_cn_pair_zmax2.txt");
    }
    {
      simprop::utils::Timer timer("timer 1");
      gnuprop::GnuProp g(std::make_unique<simprop::cosmo::Cosmology>());
      g.setRedshiftSize(1e3);
      g.setRedshiftMax(5.0);
      g.enablePhotoPion();
      g.build();
      g.evolve(0.);
      g.dump("gnuprop_spectrum_cn_full_zmax2.txt");
    }
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}