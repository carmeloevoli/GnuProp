#include "cached/PhotoPionCached.h"
#include "cached/ProtonLossesCached.h"
#include "simprop.h"

// #include "interactions/BreitWheeler.h"
// #include "interactions/KelnerAharonian2008.h"
// #include "simprop.h"
// #include "simprop/energyLosses/PairProductionLosses.h"
// #include "simprop/energyLosses/PhotoPionContinuousLosses.h"
// #include "simprop/photonFields/CmbPhotonField.h"
// #include "tables.h"

// using cmb = simprop::photonfields::CMB;

// void protonLossesTable() {
//   const auto energyAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e25 * SI::eV, 10000);
//   const auto losses_pair = simprop::losses::PairProductionLosses(std::make_unique<cmb>());
//   const auto losses_photopi =
//   simprop::losses::PhotoPionContinuousLosses(std::make_unique<cmb>()); const auto units = 1. /
//   SI::Gyr;

//   TabulatedFunction1D tabPairLosses(
//       "gnuprop_proton_pair_losses_16_25_1e4.bin",
//       [&](double energy) {
//         const auto gammap = energy / SI::protonMassC2;
//         return losses_pair.beta(simprop::proton, gammap) / units;
//       },
//       energyAxis);

//   TabulatedFunction1D tabPhotoPionLosses(
//       "gnuprop_proton_photopion_losses_16_25_1e4.bin",
//       [&](double energy) {
//         const auto gammap = energy / SI::protonMassC2;
//         return losses_photopi.beta(simprop::proton, gammap) / units;
//       },
//       energyAxis);

//   tabPairLosses.computeAndSave();
//   tabPhotoPionLosses.computeAndSave();
// }

// void nuProductionRateTable() {
//   const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 51);
//   const auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 200);
//   const auto eProtonAxis = simprop::utils::LogAxis<double>(1e19 * SI::eV, 1e22 * SI::eV, 400);
//   const auto phField = simprop::photonfields::CMB();
//   const auto nuSpec = KelnerAharonian2008::NeutrinoProductionSpectrum();
//   const auto units = 1. / SI::Gyr;

//   const std::string filename = "gnuprop_nu_production_test.bin";

//   TabulatedFunction3D tabNuProduction(
//       filename,
//       [&](double z, double x, double E_p) {
//         if (x > 1.) return 0.;

//         const auto m_pi = SI::pionMassC2;
//         const auto m_p = SI::protonMassC2;
//         const auto epsTh = (pow2(m_pi) + 2. * m_p * m_pi) / 4. / E_p;
//         const auto epsMin = std::max(epsTh, phField.getMinPhotonEnergy());
//         const auto epsMax = phField.getMaxPhotonEnergy();

//         if (epsMin > epsMax) return 0.;

//         auto integrand = [&](double lnEpsilon) {
//           const auto epsilon = std::exp(lnEpsilon);
//           const auto eta = 4. * epsilon * E_p / pow2(SI::protonMassC2);
//           auto value = phField.density(epsilon, z) * nuSpec.Phi(eta, x);
//           return epsilon * value;
//         };

//         const auto a = std::log(epsMin);
//         const auto b = std::log(epsMax);
//         const size_t N = 10000;
//         // double value = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
//         auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
//         value /= std::pow(1. + z, 3.);

//         return std::max(value / units, 0.);
//       },
//       redshifts, xAxis, eProtonAxis);

//   tabNuProduction.computeAndSave();
// }

// void gammaAbsorptionRateTable() {
//   const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 101);
//   const auto eAxis = simprop::utils::LogAxis<double>(1e13 * SI::eV, 1e22 * SI::eV, 3000);
//   const auto phField = simprop::photonfields::CMB();
//   const auto bw = PhotonPairProduction::BreitWheeler();
//   const std::string filename = "gnuprop_gamma_absorption_rate.bin";
//   const auto units = 1. / SI::Gyr;

//   TabulatedFunction2D tabGammaAbsorption(
//       filename,
//       [&](double z, double E) {
//         const auto eps_th = pow2(SI::electronMassC2) / E;
//         const auto epsMin = std::max(eps_th, phField.getMinPhotonEnergy());
//         const auto epsMax = phField.getMaxPhotonEnergy();

//         if (epsMin > epsMax) return 0.;

//         auto integrand = [&](double lnEpsilon) {
//           const auto epsilon = std::exp(lnEpsilon);
//           return phField.density(epsilon, z) / epsilon * bw.intSSigma(E * epsilon);
//         };

//         const auto a = std::log(epsMin);
//         const auto b = std::log(epsMax);
//         const size_t N = 10000;
//         auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
//         value *= SI::cLight;
//         value /= std::pow(1. + z, 3.);
//         value /= 8. * pow2(E);

//         return std::max(value / units, 0.);
//       },
//       redshifts, eAxis);

//   tabGammaAbsorption.computeAndSave();
// }

// void pairProductionRateTable() {
//   const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 51);
//   const auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 200);
//   const auto eGammaAxis = simprop::utils::LogAxis<double>(1e19 * SI::eV, 1e22 * SI::eV, 400);
//   const auto phField = simprop::photonfields::CMB();
//   const auto bw = PhotonPairProduction::BreitWheeler();
//   const auto units = 1. / SI::Gyr;

//   const std::string filename = "gnuprop_pair_production_test.bin";

//   TabulatedFunction3D tabPairProduction(
//       filename,
//       [&](double z, double x, double E_gamma) {
//         if (x > 1.) return 0.;

//         const auto E_pair = x * E_gamma;
//         const auto epsTh = 0.;  // (pow2(m_pi) + 2. * m_p * m_pi) / 4. / E_p;
//         const auto epsMin = std::max(epsTh, phField.getMinPhotonEnergy());
//         const auto epsMax = phField.getMaxPhotonEnergy();

//         if (epsMin > epsMax) return 0.;

//         auto integrand = [&](double lnEpsilon) {
//           const auto epsilon = std::exp(lnEpsilon);
//           auto value = phField.density(epsilon, z) * bw.dsigmadE(E_gamma, E_pair, epsilon);
//           return epsilon * value;
//         };

//         const auto a = std::log(epsMin);
//         const auto b = std::log(epsMax);
//         const size_t N = 10000;
//         // double value = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
//         auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
//         value /= std::pow(1. + z, 3.);

//         return std::max(value / units, 0.);
//       },
//       redshifts, xAxis, eGammaAxis);

//   tabPairProduction.computeAndSave();
// }

int main() {
  try {
    // Display startup information
    simprop::utils::startup_information();

    auto cmb = std::make_shared<simprop::photonfields::CMB>();
    auto ebl = std::make_shared<simprop::photonfields::Saldana2021PhotonField>();

    {
      simprop::utils::Timer timer("proton losses");
      const auto energyAxis = simprop::utils::LogAxis(1e16 * SI::eV, 1e25 * SI::eV, 10000);

      cache::ProtonLosses<simprop::losses::PairProductionLosses>(cmb, energyAxis,
                                                                 "gnuprop_proton_losses_pair.bin");
      cache::ProtonLosses<simprop::losses::PhotoPionContinuousLosses>(
          cmb, energyAxis, "gnuprop_proton_losses_photopion.bin");
    }
    {
      simprop::utils::Timer timer("neutrino production rate");
      const auto redshifts = simprop::utils::LinAxis<double>(0, 10, 51);
      const auto xAxis = simprop::utils::LogAxis<double>(1e-4, 1, 200);
      const auto eProtonAxis = simprop::utils::LogAxis<double>(1e19 * SI::eV, 1e22 * SI::eV, 400);

      cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(cmb, redshifts, xAxis, eProtonAxis,
                                                             "gnuprop_photopion_neutrinos_cmb.bin");
      cache::PhotoPionRate<interactions::PhotoPionNeutrinos>(ebl, redshifts, xAxis, eProtonAxis,
                                                             "gnuprop_photopion_neutrinos_ebl.bin");
    }
    // {
    //   simprop::utils::Timer timer("gamma production rate");

    //   gnuprop::PhotoPionRate<interactions::PhotoPionGammas>(cmb,
    //                                                         "gnuprop_photopion_gammas_cmb.bin");

    //   gnuprop::PhotoPionRate<interactions::PhotoPionGammas>(ebl,
    //                                                         "gnuprop_photopion_gammas_ebl.bin");
    // }
    // {
    //   simprop::utils::Timer timer("pair production rate");

    //   gnuprop::PhotoPionRate<interactions::PhotoPionPairs>(cmb,
    //   "gnuprop_photopion_pairs_cmb.bin");

    //   gnuprop::PhotoPionRate<interactions::PhotoPionPairs>(ebl,
    //   "gnuprop_photopion_pairs_ebl.bin");
    // }
    // {
    //   simprop::utils::Timer timer("gamma absorption");
    //   // gammaAbsorptionRateTable();
    // }
    // {
    //   simprop::utils::Timer timer("photon pair production");
    //   // pairProductionRateTable();
    // }
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}
