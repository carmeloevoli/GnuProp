#include "gnuprop.h"
#include "simprop.h"

int main() {
  try {
    simprop::utils::startup_information();
    {
      simprop::utils::Timer timer("timer 1");
      gnuprop::GnuProp g(std::make_unique<simprop::cosmo::Cosmology>());
      g.setRedshiftSize(1e3);
      g.setRedshiftMax(4.5);
      g.setCutoffEnergy(1e21 * SI::eV);

      g.build();
      g.evolve(0.);
      g.dump("gnuprop_spectrum_cn_db.txt");
    }

  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}