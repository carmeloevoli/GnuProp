#include "gnuprop.h"
#include "simprop.h"

// m_photoPionNu.push_back(
//     std::make_unique<gnuprop::ProductionRate>("data/gnuprop_photopion_neutrinos_cmb.bin"));
// m_photoPionNu.push_back(
//     std::make_unique<gnuprop::ProductionRate>("data/gnuprop_photopion_neutrinos_ebl.bin"));

class RunModels {
 public:
  struct ModelParameters {
    double evolutionIndex;
    double injectionSlope;
    double cutoffEnergy;
    double doCmbNeutrinos;
    double doEblNeutrinos;
    std::string outputFilename;
  };

  static std::unique_ptr<gnuprop::GnuProp> runModel(const ModelParameters& params) {
    auto g = std::make_unique<gnuprop::GnuProp>(std::make_unique<simprop::cosmo::Cosmology>());
    g->setEvolutionIndex(params.evolutionIndex);
    g->setInjectionSlope(params.injectionSlope);
    if (params.cutoffEnergy > 0) {
      g->setHeCutoff(params.cutoffEnergy);
    }
    if (params.doCmbNeutrinos) {
      g->addNeutrinoSource(
          std::make_unique<gnuprop::ProductionRate>("data/gnuprop_photopion_neutrinos_cmb.bin"));
      g->addPhotonSource(
          std::make_unique<gnuprop::ProductionRate>("data/gnuprop_photopion_gammas_cmb.bin"));
    }
    if (params.doEblNeutrinos) {
      g->addNeutrinoSource(
          std::make_unique<gnuprop::ProductionRate>("data/gnuprop_photopion_neutrinos_ebl.bin"));
      g->addPhotonSource(
          std::make_unique<gnuprop::ProductionRate>("data/gnuprop_photopion_gammas_ebl.bin"));
    }
    g->build();
    g->evolve(0.);
    g->dump(params.outputFilename);
    return g;
  }
};

int main() {
  try {
    simprop::utils::startup_information();
    {
      RunModels::ModelParameters params;
      params.evolutionIndex = 0.;
      params.injectionSlope = 2.5;
      params.cutoffEnergy = 0;
      params.doCmbNeutrinos = 1;
      params.doEblNeutrinos = 0;
      params.outputFilename = "gnuprop_spectrum_dipmodel_m0.txt";
      RunModels::runModel(params);
    }
    {
      RunModels::ModelParameters params;
      params.evolutionIndex = 3.;
      params.injectionSlope = 2.5;
      params.cutoffEnergy = 0;
      params.doCmbNeutrinos = 1;
      params.doEblNeutrinos = 0;
      params.outputFilename = "gnuprop_spectrum_dipmodel_m+3.txt";
      RunModels::runModel(params);
    }
    {
      RunModels::ModelParameters params;
      params.evolutionIndex = -3.;
      params.injectionSlope = 2.5;
      params.cutoffEnergy = 0;
      params.doCmbNeutrinos = 1;
      params.doEblNeutrinos = 0;
      params.outputFilename = "gnuprop_spectrum_dipmodel_m-3.txt";
      RunModels::runModel(params);
    }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 0.;
    //   params.injectionSlope = 2.5;
    //   params.cutoffEnergy = 7e18 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 0;
    //   params.outputFilename = "gnuprop_spectrum_auger_m0.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 3.;
    //   params.injectionSlope = 2.1;
    //   params.cutoffEnergy = 5e18 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 0;
    //   params.outputFilename = "gnuprop_spectrum_auger_m+3.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = -3.;
    //   params.injectionSlope = 2.7;
    //   params.cutoffEnergy = 5e18 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 0;
    //   params.outputFilename = "gnuprop_spectrum_auger_m-3.txt";
    //   RunModels::runModel(params);
    // }
    // EBL
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 0.;
    //   params.injectionSlope = 2.5;
    //   params.cutoffEnergy = 0;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 1;
    //   params.outputFilename = "gnuprop_spectrum_dipmodel_m0_ebl.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 3.;
    //   params.injectionSlope = 2.5;
    //   params.cutoffEnergy = 0;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 1;
    //   params.outputFilename = "gnuprop_spectrum_dipmodel_m+3_ebl.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = -3.;
    //   params.injectionSlope = 2.5;
    //   params.cutoffEnergy = 0;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 1;
    //   params.outputFilename = "gnuprop_spectrum_dipmodel_m-3_ebl.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 0.;
    //   params.injectionSlope = 2.5;
    //   params.cutoffEnergy = 7e18 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 1;
    //   params.outputFilename = "gnuprop_spectrum_auger_m0_ebl.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 3.;
    //   params.injectionSlope = 2.1;
    //   params.cutoffEnergy = 5e18 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 1;
    //   params.outputFilename = "gnuprop_spectrum_auger_m+3_ebl.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = -3.;
    //   params.injectionSlope = 2.7;
    //   params.cutoffEnergy = 5e18 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 1;
    //   params.outputFilename = "gnuprop_spectrum_auger_m-3_ebl.txt";
    //   RunModels::runModel(params);
    // }
    // {
    //   RunModels::ModelParameters params;
    //   params.evolutionIndex = 0.;
    //   params.injectionSlope = 1.3;
    //   params.cutoffEnergy = 1e20 * SI::eV;
    //   params.doCmbNeutrinos = 1;
    //   params.doEblNeutrinos = 0;
    //   params.outputFilename = "gnuprop_spectrum_uhe_m0.txt";
    //   RunModels::runModel(params);
    // }
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}