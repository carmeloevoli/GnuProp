#ifndef GNUPROP_H
#define GNUPROP_H

#include <memory>

#include "rates.h"
#include "simprop/core/cosmology.h"

namespace gnuprop {

class GnuProp {
 public:
  explicit GnuProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology);
  virtual ~GnuProp() = default;

  // parameter setters
  void setRedshiftSize(size_t N) { m_zSize = N; }
  void setRedshiftMax(double zMax) { m_zMax = zMax; }
  void setEvolutionIndex(double m) { m_evolutionIndex = m; }
  void setHeCutoff(double Ec) { m_heCutoff = Ec; }
  void setLeCutoff(double Ec) { m_leCutoff = Ec; }
  void setInjectionSlope(double slope) { m_injSlope = slope; }
  void setEnergyResolution(size_t N) { m_energySize = 13 * N; }

  // proton loss rates
  void addProtonLossRate(std::unique_ptr<gnuprop::ProtonLossRate> beta) {
    m_losses_proton.push_back(std::move(beta));
  }

  // neutrino production rates
  void addPhotoPionNeutrinoSource(std::unique_ptr<gnuprop::ProductionRate> q_nu) {
    m_sources_nu_photopion.push_back(std::move(q_nu));
  }

  // gamma production rates
  void addPhotoPionGammaSource(std::unique_ptr<gnuprop::ProductionRate> q_gamma) {
    m_sources_gamma_photopion.push_back(std::move(q_gamma));
  }
  void addInverseComptonGammaSource(std::unique_ptr<gnuprop::ProductionRate> q_gamma) {
    m_sources_gamma_ic.push_back(std::move(q_gamma));
  }
  // gamma absorption rates
  void addGammaAbsorption(std::unique_ptr<gnuprop::AbsorptionRate> k_gamma) {
    m_absorption_gamma.push_back(std::move(k_gamma));
  }
  // electron production rates
  void addPhotoPionElectronSource(std::unique_ptr<gnuprop::ProductionRate> q_electron) {
    m_sources_electron_photopion.push_back(std::move(q_electron));
  }
  void addPhotoPairElectronSource(std::unique_ptr<gnuprop::ProductionRate> q_electron) {
    m_sources_electron_photopair.push_back(std::move(q_electron));
  }
  void addInverseComptonElectronSource(std::unique_ptr<gnuprop::ProductionRate> q_electron) {
    m_sources_electrons_ic.push_back(std::move(q_electron));
  }
  // electron absorption rates
  void addElectronAbsorption(std::unique_ptr<gnuprop::AbsorptionRate> k_electron) {
    m_absorption_electron.push_back(std::move(k_electron));
  }

  // gnuprop
  void build();
  void evolve(double zObs);
  void dump(const std::string& filename) const;

 protected:
  // main parameters
  const double m_sourceComovingEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  double m_injSlope{2.5};
  double m_evolutionIndex{0.0};
  double m_heCutoff{1e23 * SI::eV};
  double m_leCutoff{1e17 * SI::eV};
  double m_zMax{5.0};
  const double m_energyMin{1e10 * SI::eV};
  const double m_energyMax{1e23 * SI::eV};
  size_t m_energySize{13 * 32};
  size_t m_zSize{100000};

 protected:
  // cosmology
  std::unique_ptr<simprop::cosmo::Cosmology> m_cosmology;

  // energy axis
  std::vector<double> m_eAxis;

  // protons
  std::vector<std::unique_ptr<gnuprop::ProtonLossRate>> m_losses_proton;
  std::vector<double> m_n_proton;
  std::vector<double> m_q_proton;
  std::vector<double> m_beta_proton;

  // neutrinos
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_nu_photopion;
  std::vector<double> m_q_nu;
  std::vector<double> m_n_nu;

  // photons
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_gamma_photopion;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_gamma_ic;
  std::vector<std::unique_ptr<gnuprop::AbsorptionRate>> m_absorption_gamma;
  std::vector<double> m_q_gamma;
  std::vector<double> m_n_gamma;
  std::vector<double> m_k_gamma;

  // electrons
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_electron_photopion;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_electron_photopair;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_electrons_ic;
  std::vector<std::unique_ptr<gnuprop::AbsorptionRate>> m_absorption_electron;
  std::vector<double> m_q_electron;
  std::vector<double> m_n_electron;
  std::vector<double> m_k_electron;

  // CN method temporary vectors
  std::vector<double> knownTerm;
  std::vector<double> diagonal;
  std::vector<double> upperDiagonal;
  std::vector<double> lowerDiagonal;

 protected:
  // source functions
  void evolveProtonEmissivity(double z);
  void evolveProtonLosses(double z);
  void evolveNuEmissivity(double z);
  void evolveGammaEmissivity(double z);
  void evolveGammaAbsorption(double z);
  void evolveElectronEmissivity(double z);
  void evolveElectronAbsorption(double z);
};

}  // namespace gnuprop

#endif  // GAMMAPROP_H