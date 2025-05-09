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
  // template <typename T>
  //  static std::unique_ptr<cache::ProtonLossesCached> protonLosses() {
  //    auto losses = std::make_unique<cache::ProtonLossesCached>();
  //    losses->buildEnergyAxis(1e16 * SI::eV, 1e24 * SI::eV, 8 * 64);
  //    losses->builLosses(std::make_unique<T>(cmb));
  //    return losses;
  //  }

  // template <typename T>
  // static std::unique_ptr<T> absorptionCached(
  //     const std::shared_ptr<simprop::photonfields::PhotonField>& photonField) {
  //   auto absorption = std::make_unique<T>();
  //   absorption->buildPrimaryEnergyAxis(1e10 * SI::eV, 1e24 * SI::eV, 14 * 32);
  //   absorption->buildRedshiftAxis(0., 5., 101);
  //   absorption->buildPhotonField(photonField);
  //   return absorption;
  // }

  template <typename T>
  static std::unique_ptr<cache::PhotoPionCached<T>> photoPionCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField) {
    auto photoPion = std::make_unique<cache::PhotoPionCached<T>>();
    photoPion->buildEnergyAxis(1e15 * SI::eV, 1e24 * SI::eV, 9 * 16);
    photoPion->buildXAxis(1e-6, 1, 6 * 64);
    photoPion->buildRedshiftAxis(0., 5., 101);
    photoPion->buildPhotonField(photonField);
    return photoPion;
  }

  static std::unique_ptr<cache::PhotoPairCached> photoPairCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField) {
    auto photoPair = std::make_unique<cache::PhotoPairCached>();
    photoPair->buildEnergyAxis(1e10 * SI::eV, 1e24 * SI::eV, 14 * 32);
    photoPair->buildXAxis(1e-5, 1, 5 * 100);
    photoPair->buildRedshiftAxis(0., 5., 101);
    photoPair->buildPhotonField(photonField);
    return photoPair;
  }

  static std::unique_ptr<cache::InverseComptonCached> inverseComptonCached(
      const std::shared_ptr<simprop::photonfields::PhotonField>& photonField,
      bool doGammas = true) {
    auto ic = std::make_unique<cache::InverseComptonCached>();
    ic->buildEnergyAxis(1e10 * SI::eV, 1e24 * SI::eV, 14 * 32);
    ic->buildXAxis(1e-5, 1, 5 * 100);
    ic->buildRedshiftAxis(0., 5., 101);
    ic->buildPhotonField(photonField);
    ic->setDoGammas(doGammas);
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
    // gammaAbsCmb->run("gnuprop_absorption_gammas_cmb.bin", 1e-3);
    // auto pairAbsCmb = GnuPropTables::absorptionCached<cache::PairAbsorptionCached>(cmb);
    // pairAbsCmb->run("gnuprop_absorption_pairs_cmb.bin", 1e-3);
    // auto gammaAbsEbl = GnuPropTables::absorptionCached<cache::GammaAbsorptionCached>(ebl);
    // gammaAbsEbl->run("gnuprop_absorption_gammas_ebl.bin", 1e-2);
    // auto pairAbsEbl = GnuPropTables::absorptionCached<cache::PairAbsorptionCached>(ebl);
    // pairAbsEbl->run("gnuprop_absorption_pairs_ebl.bin", 1e-2);

    // photo-pion rates
    auto ppNuCmb = GnuPropTables::photoPionCached<interactions::PhotoPionNeutrinos>(cmb);
    ppNuCmb->run("gnuprop_photopion_neutrinos_cmb.bin", 1e-3);
    auto ppNuEbl = GnuPropTables::photoPionCached<interactions::PhotoPionNeutrinos>(ebl);
    ppNuEbl->run("gnuprop_photopion_neutrinos_ebl.bin", 1e-2);
    auto ppPairCmb = GnuPropTables::photoPionCached<interactions::PhotoPionPairs>(cmb);
    ppPairCmb->run("gnuprop_photopion_pairs_cmb.bin", 1e-3);
    auto ppPairEbl = GnuPropTables::photoPionCached<interactions::PhotoPionPairs>(ebl);
    ppPairEbl->run("gnuprop_photopion_pairs_ebl.bin", 1e-2);
    auto ppGammaCmb = GnuPropTables::photoPionCached<interactions::PhotoPionGammas>(cmb);
    ppGammaCmb->run("gnuprop_photopion_gammas_cmb.bin", 1e-3);
    auto ppGammaEbl = GnuPropTables::photoPionCached<interactions::PhotoPionGammas>(ebl);
    ppGammaEbl->run("gnuprop_photopion_gammas_ebl.bin", 1e-2);

    // photo-pair rates
    // auto ppCmb = GnuPropTables::photoPairCached(cmb);
    // ppCmb->run("gnuprop_photopair_pairs_cmb.bin", 1e-3);
    // auto ppEbl = GnuPropTables::photoPairCached(ebl);
    // ppEbl->run("gnuprop_photopair_pairs_ebl.bin", 1e-2);

    // IC rates
    // auto icCmb = GnuPropTables::inverseComptonCached(cmb);
    // icCmb->run("gnuprop_inversecompton_gammas_cmb.bin", 1e-3);
    // auto icEbl = GnuPropTables::inverseComptonCached(ebl);
    //  icEbl->run("gnuprop_inversecompton_gammas_ebl.bin", 1e-2);
    //  auto icPairCmb = GnuPropTables::inverseComptonCached(cmb, false);
    //  icPairCmb->run("gnuprop_inversecompton_pairs_cmb.bin", 1e-3);
    //  auto icPairEbl = GnuPropTables::inverseComptonCached(ebl, false);
    //  icPairEbl->run("gnuprop_inversecompton_pairs_ebl.bin", 1e-2);
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}