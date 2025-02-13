#include "cached/GammaAbsorptionCached.h"
#include "cached/InverseComptonCached.h"
#include "cached/PairAbsorptionCached.h"
#include "cached/PhotoPairCached.h"
#include "cached/PhotoPionCached.h"
#include "cached/ProtonLossesCached.h"
#include "simprop.h"

const auto cmb = std::make_shared<simprop::photonfields::CMB>();
const auto ebl = std::make_shared<simprop::photonfields::Saldana2021PhotonField>();
const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 51);

void compute_proton_losses() {
  simprop::utils::Timer timer("proton losses");
  const auto energyAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e23 * SI::eV, 1000);

  cache::ProtonLosses<simprop::losses::PairProductionLosses>(cmb, energyAxis,
                                                             "gnuprop_proton_losses_pair.bin");
  cache::ProtonLosses<simprop::losses::PhotoPionContinuousLosses>(
      cmb, energyAxis, "gnuprop_proton_losses_photopion.bin");
}

void compute_photopair_absorption_rates() {
  simprop::utils::Timer timer("photo-pair absorption rates");

  const auto energyAxis = simprop::utils::LogAxis<double>(1e10 * SI::eV, 1e23 * SI::eV, 1000);

  cache::GammaAbsorptionRate(cmb, redshifts, energyAxis, "gnuprop_absorption_gammas_cmb.bin");
  cache::GammaAbsorptionRate(ebl, redshifts, energyAxis, "gnuprop_absorption_gammas_ebl.bin");
}

void compute_ic_absorption_rates() {
  simprop::utils::Timer timer("IC absorption rates");

  const auto energyAxis = simprop::utils::LogAxis<double>(1e10 * SI::eV, 1e23 * SI::eV, 1000);

  cache::PairAbsorptionRate(cmb, redshifts, energyAxis, "gnuprop_absorption_pairs_cmb.bin");
  cache::PairAbsorptionRate(ebl, redshifts, energyAxis, "gnuprop_absorption_pairs_ebl.bin");
}

void compute_photo_pion_rates() {
  simprop::utils::Timer timer("photo-pion production rate");

  const auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 100);
  const auto eProtonAxis = simprop::utils::LogAxis<double>(1e17 * SI::eV, 1e23 * SI::eV, 1000);

  cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(ebl, redshifts, xAxis, eProtonAxis,
                                                         "gnuprop_photopion_neutrinos_ebl.bin");
  cache::PhotoPionRate<interactions::PhotoPionGammas>(ebl, redshifts, xAxis, eProtonAxis,
                                                      "gnuprop_photopion_gammas_ebl.bin");
  cache::PhotoPionRate<interactions::PhotoPionPairs>(ebl, redshifts, xAxis, eProtonAxis,
                                                     "gnuprop_photopion_pairs_ebl.bin");
  cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(cmb, redshifts, xAxis, eProtonAxis,
                                                         "gnuprop_photopion_neutrinos_cmb.bin");
  cache::PhotoPionRate<interactions::PhotoPionGammas>(cmb, redshifts, xAxis, eProtonAxis,
                                                      "gnuprop_photopion_gammas_cmb.bin");
  cache::PhotoPionRate<interactions::PhotoPionPairs>(cmb, redshifts, xAxis, eProtonAxis,
                                                     "gnuprop_photopion_pairs_cmb.bin");
}

void compute_photo_pair_rates() {
  simprop::utils::Timer timer("photo-pair production rate");

  const auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 100);
  const auto eGammaAxis = simprop::utils::LogAxis<double>(1e17 * SI::eV, 1e23 * SI::eV, 1000);

  cache::PhotoPairRate(cmb, redshifts, xAxis, eGammaAxis, "gnuprop_photopair_pairs_cmb.bin");
  cache::PhotoPairRate(ebl, redshifts, xAxis, eGammaAxis, "gnuprop_photopair_pairs_ebl.bin");
}

void compute_IC_rates() {
  simprop::utils::Timer timer("IC production rate");

  const auto xAxis = simprop::utils::LogAxis<double>(1e-5, 1, 100);
  const auto eElectronAxis = simprop::utils::LogAxis<double>(1e17 * SI::eV, 1e23 * SI::eV, 1000);

  cache::InverseComptonRate(cmb, redshifts, xAxis, eElectronAxis, "gnuprop_IC_gammas_cmb.bin");
  cache::InverseComptonRate(ebl, redshifts, xAxis, eElectronAxis, "gnuprop_IC_gammas_ebl.bin");
  cache::InverseComptonRate(cmb, redshifts, xAxis, eElectronAxis, "gnuprop_IC_pairs_cmb.bin",
                            false);
  cache::InverseComptonRate(ebl, redshifts, xAxis, eElectronAxis, "gnuprop_IC_pairs_ebl.bin",
                            false);
}

int main() {
  try {
    // Display startup information
    simprop::utils::startup_information();
    // + proton losses
    // compute_proton_losses();
    // + absorption rates
    // compute_photopair_absorption_rates();
    // compute_ic_absorption_rates();
    // + production rates
    // compute_photo_pion_rates();
    // compute_photo_pair_rates();
    // compute_IC_rates();
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}
