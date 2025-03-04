#include "cached/GammaAbsorptionCached.h"
#include "cached/InverseComptonCached.h"
#include "cached/PairAbsorptionCached.h"
#include "cached/PhotoPairCached.h"
#include "cached/PhotoPionCached.h"
#include "cached/ProtonLossesCached.h"
#include "simprop.h"

const auto cmb = std::make_shared<simprop::photonfields::CMB>();
const auto ebl = std::make_shared<simprop::photonfields::Saldana2021PhotonField>();

class GnuPropTables {
 public:
  template <typename T>
  static std::unique_ptr<cache::ProtonLossesCached> protonLosses() {
    auto losses = std::make_unique<cache::ProtonLossesCached>();
    losses->buildEnergyAxis(1e16 * SI::eV, 1e24 * SI::eV, 8 * 128);
    losses->builLosses(std::make_unique<T>(cmb));
    return losses;
  }

  template <typename T>
  static std::unique_ptr<T> absorptionCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField) {
    auto absorption = std::make_unique<T>();
    absorption->buildEnergyAxis(1e10 * SI::eV, 1e24 * SI::eV, 14 * 128);
    absorption->buildRedshiftAxis(0., 5., 101);
    absorption->buildPhotonField(photonField);
    return absorption;
  }

  template <typename T>
  static std::unique_ptr<cache::PhotoPionCached<T>> photoPionCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField,
      double precision = 1e-4) {
    auto photoPion = std::make_unique<cache::PhotoPionCached<T>>();
    photoPion->buildEnergyAxis(1e16 * SI::eV, 1e24 * SI::eV, 8 * 64);
    photoPion->buildXAxis(1e-5, 1., 5 * 64);
    photoPion->buildRedshiftAxis(0., 5., 101);
    photoPion->buildPhotonField(photonField);
    photoPion->setPrecision(precision);
    return photoPion;
  }

  static std::unique_ptr<cache::PhotoPairCached> photoPairCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField,
      double precision = 1e-4) {
    auto photoPair = std::make_unique<cache::PhotoPairCached>();
    photoPair->buildEnergyAxis(1e12 * SI::eV, 1e24 * SI::eV, 12 * 64);
    photoPair->buildXAxis(1e-5, 1., 5 * 64);
    photoPair->buildRedshiftAxis(0., 5., 101);
    photoPair->buildPhotonField(photonField);
    photoPair->setPrecision(precision);
    return photoPair;
  }

  static std::unique_ptr<cache::InverseComptonCached> inverseComptonCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField,
      bool doGammas = true) {
    auto ic = std::make_unique<cache::InverseComptonCached>();
    ic->buildEnergyAxis(1e12 * SI::eV, 1e24 * SI::eV, 12 * 64);
    ic->buildXAxis(1e-5, 1., 5 * 64);
    ic->buildRedshiftAxis(0., 5., 101);
    ic->setDoGammas(doGammas);
    ic->buildPhotonField(photonField);
    return ic;
  }
};

int main() {
  try {
    // Display startup information
    simprop::utils::startup_information();

    // proton losses
    // auto pairLosses = GnuPropTables::protonLosses<simprop::losses::PairProductionLosses>();
    // pairLosses->run("gnuprop_proton_losses_pair.bin");
    // auto pionLosses = GnuPropTables::protonLosses<simprop::losses::PhotoPionContinuousLosses>();
    // pionLosses->run("gnuprop_proton_losses_photopion.bin");

    // absorption rates
    // auto gammaAbsCmb = GnuPropTables::absorptionCached<cache::GammaAbsorptionCached>(cmb);
    // gammaAbsCmb->run("gnuprop_absorption_gammas_cmb.bin");
    // auto gammaAbsEbl = GnuPropTables::absorptionCached<cache::GammaAbsorptionCached>(ebl);
    // gammaAbsEbl->run("gnuprop_absorption_gammas_ebl.bin");
    // auto pairAbsCmb = GnuPropTables::absorptionCached<cache::PairAbsorptionCached>(cmb);
    // pairAbsCmb->run("gnuprop_absorption_pairs_cmb.bin");
    // auto pairAbsEbl = GnuPropTables::absorptionCached<cache::PairAbsorptionCached>(ebl);
    // pairAbsEbl->run("gnuprop_absorption_pairs_ebl.bin");

    // photo-pion rates
    // auto ppNuCmb = GnuPropTables::photoPionCached<interactions::PhotoPionNeutrinos>(cmb);
    // ppNuCmb->run("gnuprop_photopion_neutrinos_cmb.bin");
    // auto ppGammaCmb = GnuPropTables::photoPionCached<interactions::PhotoPionGammas>(cmb);
    // ppGammaCmb->run("gnuprop_photopion_gammas_cmb.bin");
    // auto ppPairCmb = GnuPropTables::photoPionCached<interactions::PhotoPionPairs>(cmb);
    // ppPairCmb->run("gnuprop_photopion_pairs_cmb.bin");
    // auto ppNuEbl = GnuPropTables::photoPionCached<interactions::PhotoPionNeutrinos>(ebl);
    // ppNuEbl->run("gnuprop_photopion_neutrinos_ebl.bin");
    // auto ppGammaEbl = GnuPropTables::photoPionCached<interactions::PhotoPionGammas>(ebl);
    // ppGammaEbl->run("gnuprop_photopion_gammas_ebl.bin");
    // auto ppPairEbl = GnuPropTables::photoPionCached<interactions::PhotoPionPairs>(ebl);
    // ppPairEbl->run("gnuprop_photopion_pairs_ebl.bin");

    // photo-pair rates
    // auto ppCmb = GnuPropTables::photoPairCached(cmb);
    // ppCmb->run("gnuprop_photopair_pairs_cmb.bin");
    // auto ppEbl = GnuPropTables::photoPairCached(ebl, 1e-2);
    // ppEbl->run("gnuprop_photopair_pairs_ebl.bin");

    // IC rates
    auto icCmb = GnuPropTables::inverseComptonCached(cmb);
    icCmb->run("gnuprop_inversecompton_gammas_cmb.bin");
    auto icEbl = GnuPropTables::inverseComptonCached(ebl);
    icEbl->run("gnuprop_inversecompton_gammas_ebl.bin");
    auto icPairCmb = GnuPropTables::inverseComptonCached(cmb, false);
    icPairCmb->run("gnuprop_inversecompton_pairs_cmb.bin");
    auto icPairEbl = GnuPropTables::inverseComptonCached(ebl, false);
    icPairEbl->run("gnuprop_inversecompton_pairs_ebl.bin");
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}

// void compute_photo_pair_rates() {
//   simprop::utils::Timer timer("photo-pair production rate");

//   cache::PhotoPairRate(ebl, redshifts, xAxis, energyAxis, "gnuprop_photopair_pairs_ebl.bin");
//   cache::PhotoPairRate(cmb, redshifts, xAxis, energyAxis, "gnuprop_photopair_pairs_cmb.bin");
// }

// compute_photo_pion_rates();
//   compute_photo_pair_rates();
//   compute_IC_rates();

// static std::unique_ptr<cache::PairAbsorptionCached> pairAbsorptionCmb() {
//   auto ic = std::make_unique<cache::PairAbsorptionCached>();
//   ic->buildEnergyAxis(1e15 * SI::eV, 1e22 * SI::eV, 1000);
//   ic->buildRedshiftAxis(0., 6., 6);
//   ic->buildPhotonField(cmb);
//   return ic;
// }

// static std::unique_ptr<cache::PairAbsorptionCached> pairAbsorptionEbl() {
//   auto ic = std::make_unique<cache::PairAbsorptionCached>();
//   ic->buildEnergyAxis(1e15 * SI::eV, 1e22 * SI::eV, 1000);
//   ic->buildRedshiftAxis(0., 6., 6);
//   ic->buildPhotonField(ebl);
//   return ic;
// }

// void compute_proton_losses() {
//   simprop::utils::Timer timer("proton losses");

//   const auto eProtonAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e24 * SI::eV, 8 * 128);

//   cache::ProtonLosses<simprop::losses::PairProductionLosses>(cmb, eProtonAxis,
//                                                              "gnuprop_proton_losses_pair.bin");
//   cache::ProtonLosses<simprop::losses::PhotoPionContinuousLosses>(
//       cmb, eProtonAxis, "gnuprop_proton_losses_photopion.bin");
// }

// void compute_photopair_absorption_rates() {
//   simprop::utils::Timer timer("photo-pair absorption rates");

//   cache::GammaAbsorptionRate(cmb, redshifts, energyAxis, "gnuprop_absorption_gammas_cmb.bin");
//   cache::GammaAbsorptionRate(ebl, redshifts, energyAxis, "gnuprop_absorption_gammas_ebl.bin");
// }

// void compute_ic_absorption_rates() {
//   simprop::utils::Timer timer("IC absorption rates");

//   cache::PairAbsorptionRate(cmb, redshifts, energyAxis, "gnuprop_absorption_pairs_cmb.bin");
//   cache::PairAbsorptionRate(ebl, redshifts, energyAxis, "gnuprop_absorption_pairs_ebl.bin");
// }

// void compute_photo_pion_rates() {
//   simprop::utils::Timer timer("photo-pion production rate");

//   cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(ebl, redshifts, xAxis, energyAxis,
//                                                          "gnuprop_photopion_neutrinos_ebl.bin");
//   cache::PhotoPionRate<interactions::PhotoPionGammas>(ebl, redshifts, xAxis, energyAxis,
//                                                       "gnuprop_photopion_gammas_ebl.bin");
//   cache::PhotoPionRate<interactions::PhotoPionPairs>(ebl, redshifts, xAxis, energyAxis,
//                                                      "gnuprop_photopion_pairs_ebl.bin");
//   cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(cmb, redshifts, xAxis, energyAxis,
//                                                          "gnuprop_photopion_neutrinos_cmb.bin");
//   cache::PhotoPionRate<interactions::PhotoPionGammas>(cmb, redshifts, xAxis, energyAxis,
//                                                       "gnuprop_photopion_gammas_cmb.bin");
//   cache::PhotoPionRate<interactions::PhotoPionPairs>(cmb, redshifts, xAxis, energyAxis,
//                                                      "gnuprop_photopion_pairs_cmb.bin");
// }

// void compute_photo_pair_rates() {
//   simprop::utils::Timer timer("photo-pair production rate");

//   cache::PhotoPairRate(ebl, redshifts, xAxis, energyAxis, "gnuprop_photopair_pairs_ebl.bin");
//   cache::PhotoPairRate(cmb, redshifts, xAxis, energyAxis, "gnuprop_photopair_pairs_cmb.bin");
// }

// void compute_IC_rates() {
//   simprop::utils::Timer timer("IC production rate");

//   cache::InverseComptonRate(cmb, redshifts, xAxis, energyAxis, "gnuprop_IC_gammas_cmb.bin");
//   cache::InverseComptonRate(ebl, redshifts, xAxis, energyAxis, "gnuprop_IC_gammas_ebl.bin");
//   cache::InverseComptonRate(cmb, redshifts, xAxis, energyAxis, "gnuprop_IC_pairs_cmb.bin",
//   false); cache::InverseComptonRate(ebl, redshifts, xAxis, energyAxis,
//   "gnuprop_IC_pairs_ebl.bin", false);
// }