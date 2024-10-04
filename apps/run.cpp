#include "gammaprop.h"
#include "simprop.h"

int main() {
  try {
    simprop::utils::startup_information();
    simprop::utils::Timer timer("main timer");
    {
      gammaprop::GammaProp g(std::make_unique<simprop::cosmo::Cosmology>());
      g.evolve(0.);
      g.dump("GammaProp_spectrum_full_z0.txt");
    }

    // {
    //   gammaprop::GammaProp g(std::make_unique<simprop::cosmo::Cosmology>());
    //   g.evolve(1.);
    //   g.dump("GammaProp_spectrum_full_z1.txt");
    // }

    // {
    //   gammaprop::GammaProp g(std::make_unique<simprop::cosmo::Cosmology>());
    //   g.evolve(2.);
    //   g.dump("GammaProp_spectrum_full_z2.txt");
    // }

  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}