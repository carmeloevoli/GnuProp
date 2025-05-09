#include "gnuprop.h"
#include "simprop.h"

class RunModels {
 public:
  struct ModelParameters {
    double evolutionIndex = 0.;
    double zMax = 5.;
    double redshiftSize = 10000;  // 100000
    double injectionSlope = 2.5;
    double cutoffEnergy = 0.;
    bool doCmbNeutrinos = false;
    bool doEblNeutrinos = false;
    std::string outputFilename;
  };

  static std::unique_ptr<gnuprop::GnuProp> runModel(const ModelParameters& params) {
    simprop::utils::Timer timer("Time for running model " + params.outputFilename);

    auto g = std::make_unique<gnuprop::GnuProp>(std::make_unique<simprop::cosmo::Cosmology>());
    g->setRedshiftMax(params.zMax);
    g->setRedshiftSize(params.redshiftSize);
    g->setEvolutionIndex(params.evolutionIndex);
    g->setInjectionSlope(params.injectionSlope);
    if (params.cutoffEnergy > 0) {
      g->setHeCutoff(params.cutoffEnergy);
    }

    g->addProtonLossRate(
        std::make_unique<gnuprop::ProtonLossRate>("tables/gnuprop_proton_losses_pair.bin"));
    g->addProtonLossRate(
        std::make_unique<gnuprop::ProtonLossRate>("tables/gnuprop_proton_losses_photopion.bin"));

    if (params.doCmbNeutrinos) {
      g->addPhotoPionNeutrinoSource(
          std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_neutrinos_cmb.bin"));
      g->addPhotoPionGammaSource(
          std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_gammas_cmb.bin"));
      g->addPhotoPionElectronSource(
          std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_pairs_cmb.bin"));

      // g->addGammaAbsorption(
      //     std::make_unique<gnuprop::AbsorptionRate>("tables/gnuprop_absorption_gammas_cmb.bin"));
      // g->addPhotoPairElectronSource(
      //     std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopair_pairs_cmb.bin"));

      // g->addElectronAbsorption(
      //    std::make_unique<gnuprop::AbsorptionRate>("tables/gnuprop_absorption_pairs_cmb.bin"));
      // g->addInverseComptonGammaSource(std::make_unique<gnuprop::ProductionRate>(
      //     "tables/gnuprop_inversecompton_gammas_cmb.bin"));
      // g->addInverseComptonElectronSource(
      //     std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_inversecompton_pairs_cmb.bin"));
    }

    if (params.doEblNeutrinos) {
      g->addPhotoPionNeutrinoSource(
          std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_neutrinos_ebl.bin"));
      g->addPhotoPionGammaSource(
          std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_gammas_ebl.bin"));
      g->addPhotoPionElectronSource(
          std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopion_pairs_ebl.bin"));

      // g->addGammaAbsorption(
      //     std::make_unique<gnuprop::AbsorptionRate>("tables/gnuprop_absorption_gammas_ebl.bin"));

      // g->addInverseComptonGammaSource(std::make_unique<gnuprop::ProductionRate>(
      //     "tables/gnuprop_inversecompton_gammas_ebl.bin"));

      // g->addPhotoPairElectronSource(
      //     std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_photopair_pairs_ebl.bin"));
      // g->addInverseComptonElectronSource(
      //     std::make_unique<gnuprop::ProductionRate>("tables/gnuprop_inversecompton_pairs_ebl.bin"));
    }
    g->build();
    g->evolve(0.);
    g->dump(params.outputFilename);
    return g;
  }
};

void run_dip() {
  const double MAXCUTOFF = 1e21 * SI::eV;

  // Dip Model
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 0.;
    params.injectionSlope = 2.6;
    params.cutoffEnergy = MAXCUTOFF;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_dipmodel_m0_cmb.txt";
    RunModels::runModel(params);
    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_dipmodel_m0_ebl.txt";
    RunModels::runModel(params);
  }
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 3.;
    params.injectionSlope = 2.5;
    params.cutoffEnergy = MAXCUTOFF;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_dipmodel_m3_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_dipmodel_m3_ebl.txt";
    RunModels::runModel(params);
  }
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 5.;
    params.injectionSlope = 2.4;
    params.cutoffEnergy = MAXCUTOFF;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_dipmodel_m5_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_dipmodel_m5_ebl.txt";
    RunModels::runModel(params);
  }
}

void run_auger() {
  // Auger Proton Model
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 0.;
    params.injectionSlope = 2.5;
    params.cutoffEnergy = 6e18 * SI::eV;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_auger_m0_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_auger_m0_ebl.txt";
    RunModels::runModel(params);
  }
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 3.;
    params.injectionSlope = 2.01;
    params.cutoffEnergy = 5e18 * SI::eV;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_auger_m3_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_auger_m3_ebl.txt";
    RunModels::runModel(params);
  }
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 5.;
    params.injectionSlope = 1.4;
    params.cutoffEnergy = 4e18 * SI::eV;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_auger_m5_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_auger_m5_ebl.txt";
    RunModels::runModel(params);
  }
}

void run_uhe() {
  const double MAXCUTOFF = 1e20 * SI::eV;

  // UHE Model
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 0.;
    params.injectionSlope = 1.3;
    params.cutoffEnergy = MAXCUTOFF;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_uhe_m0_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_uhe_m0_ebl.txt";
    RunModels::runModel(params);
  }
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 3.;
    params.injectionSlope = 1.0;
    params.cutoffEnergy = MAXCUTOFF;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_uhe_m3_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_uhe_m3_ebl.txt";
    RunModels::runModel(params);
  }
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 5.;
    params.injectionSlope = 0.7;
    params.cutoffEnergy = MAXCUTOFF;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_uhe_m5_cmb.txt";
    RunModels::runModel(params);

    params.doEblNeutrinos = true;
    params.outputFilename = "gnuprop_spectrum_uhe_m5_ebl.txt";
    RunModels::runModel(params);
  }
}

void run_test() {
  {
    RunModels::ModelParameters params;
    params.evolutionIndex = 3.;
    params.injectionSlope = 2.5;
    params.zMax = 5.;
    params.cutoffEnergy = 1e23 * SI::eV;
    params.doCmbNeutrinos = true;
    params.doEblNeutrinos = false;
    params.outputFilename = "gnuprop_spectrum_test_m3_cmb.txt";
    RunModels::runModel(params);
  }
}

int main() {
  try {
    simprop::utils::startup_information();

    run_uhe();

  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}
