#include "cached/GammaAbsorptionCached.h"
#include "cached/PairAbsorptionCached.h"
#include "cached/PhotoPionCached.h"
#include "cached/ProtonLossesCached.h"
#include "simprop.h"

void protonlosses() {
  simprop::utils::Timer timer("proton losses");

  auto cmb = std::make_shared<simprop::photonfields::CMB>();
  const auto energyAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e24 * SI::eV, 400);

  cache::ProtonLosses<simprop::losses::PairProductionLosses>(cmb, energyAxis,
                                                             "gnuprop_proton_losses_pair.bin");
  cache::ProtonLosses<simprop::losses::PhotoPionContinuousLosses>(
      cmb, energyAxis, "gnuprop_proton_losses_photopion.bin");
}

void photopionrates() {
  simprop::utils::Timer timer("photo-pion production rate");
  const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 21);
  const auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 60);
  const auto eProtonAxis = simprop::utils::LogAxis<double>(1e17 * SI::eV, 1e23 * SI::eV, 200);

  auto cmb = std::make_shared<simprop::photonfields::CMB>();
  auto ebl = std::make_shared<simprop::photonfields::Saldana2021PhotonField>();

  cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(cmb, redshifts, xAxis, eProtonAxis,
                                                         "gnuprop_photopion_neutrinos_cmb.bin");
  cache::PhotoPionRate<interactions::PhotoPionGammas>(cmb, redshifts, xAxis, eProtonAxis,
                                                      "gnuprop_photopion_gammas_cmb.bin");
  cache::PhotoPionRate<interactions::PhotoPionPairs>(cmb, redshifts, xAxis, eProtonAxis,
                                                     "gnuprop_photopion_pairs_cmb.bin");
  cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(ebl, redshifts, xAxis, eProtonAxis,
                                                         "gnuprop_photopion_neutrinos_ebl.bin");
  cache::PhotoPionRate<interactions::PhotoPionGammas>(ebl, redshifts, xAxis, eProtonAxis,
                                                      "gnuprop_photopion_gammas_ebl.bin");
  cache::PhotoPionRate<interactions::PhotoPionPairs>(ebl, redshifts, xAxis, eProtonAxis,
                                                     "gnuprop_photopion_pairs_ebl.bin");
}

void absorptionrates() {
  simprop::utils::Timer timer("absorption rates");

  auto cmb = std::make_shared<simprop::photonfields::CMB>();
  auto ebl = std::make_shared<simprop::photonfields::Saldana2021PhotonField>();

  const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 101);
  const auto energyAxis = simprop::utils::LogAxis<double>(1e12 * SI::eV, 1e22 * SI::eV, 1000);

  cache::GammaAbsorptionRate(cmb, redshifts, energyAxis, "gnuprop_absorption_gammas_cmb.bin");
  cache::PairAbsorptionRate(cmb, redshifts, energyAxis, "gnuprop_absorption_pairs_cmb.bin");
  cache::GammaAbsorptionRate(ebl, redshifts, energyAxis, "gnuprop_absorption_gammas_ebl.bin");
  cache::PairAbsorptionRate(ebl, redshifts, energyAxis, "gnuprop_absorption_pairs_ebl.bin");
}

int main() {
  try {
    // Display startup information
    simprop::utils::startup_information();
    // protonlosses();
    // photopionrates();
    absorptionrates();
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}
