#ifndef GNUPROP_H
#define GNUPROP_H

#include <memory>
#include <stdexcept>

#include "rates.h"
#include "simprop/core/cosmology.h"

namespace gnuprop {

template <typename T>
double compute_integral(const std::vector<T> &x, const std::vector<T> &f) {
  if (x.size() != f.size()) {
    throw std::invalid_argument("x and f must have the same size.");
  }
  if (x.size() < 2) {
    throw std::invalid_argument("At least two points are required.");
  }

  T integral = 0.0;
  for (std::size_t i = 0; i < x.size() - 1; ++i) {
    const T dx = x[i + 1] - x[i];
    integral += 0.5 * (f[i] + f[i + 1]) * dx;
  }

  return integral;
}

template <typename T>
double compute_energy_integral(const std::vector<T> &x, const std::vector<T> &f) {
  if (x.size() != f.size()) {
    throw std::invalid_argument("x and f must have the same size.");
  }
  if (x.size() < 2) {
    throw std::invalid_argument("At least two points are required.");
  }

  T integral = 0.0;
  for (std::size_t i = 0; i < x.size() - 1; ++i) {
    const T dx = x[i + 1] - x[i];
    integral += 0.5 * (x[i] * f[i] + x[i + 1] * f[i + 1]) * dx;
  }

  return integral;
}

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
  void dump(const std::string &filename) const;

 public:
  // main parameters
  const double m_sourceComovingEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  const double E_0 = 1e18 * SI::eV;
  double m_sourceNormalization{0.0};
  double m_injSlope{2.5};
  double m_evolutionIndex{0.0};
  double m_heCutoff{1e23 * SI::eV};
  double m_leCutoff{1e17 * SI::eV};
  double m_zMax{5.0};
  double m_dz{0.0};
  const double m_energyMin{1e10 * SI::eV};
  const double m_energyMax{1e23 * SI::eV};
  size_t m_energySize{13 * 32};
  size_t m_zSize{10000};

 public:
  // cosmology
  std::unique_ptr<simprop::cosmo::Cosmology> m_cosmology;

  // source functions
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_nu_photopion;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_gamma_photopion;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_gamma_ic;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_electron_photopion;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_electron_photopair;
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_sources_electrons_ic;

  // losses functions
  std::vector<std::unique_ptr<gnuprop::ProtonLossRate>> m_losses_proton;
  std::vector<std::unique_ptr<gnuprop::AbsorptionRate>> m_absorption_electron;
  std::vector<std::unique_ptr<gnuprop::AbsorptionRate>> m_absorption_gamma;

  // energy axis
  std::vector<double> m_eAxis;

  // protons
  std::vector<double> m_q_proton;
  std::vector<double> m_beta_proton;
  std::vector<double> m_f_proton;

  // neutrinos
  std::vector<double> m_q_nu;
  std::vector<double> m_f_nu;

  // photons
  std::vector<double> m_q_gamma;
  std::vector<double> m_k_gamma;
  std::vector<double> m_f_gamma;

  // electrons
  std::vector<double> m_q_electron;
  std::vector<double> m_k_electron;
  std::vector<double> m_f_electron;

  // CN method temporary vectors
  std::vector<double> knownTerm;
  std::vector<double> diagonal;
  std::vector<double> upperDiagonal;
  std::vector<double> lowerDiagonal;

 public:
  // source functions
  void computeProtonEmissivity(double z);
  void computeProtonLosses(double z);
  void computeNuEmissivity(double z);
  void computeGammaEmissivity(double z);
  void computeElectronEmissivity(double z);
  void computeGammaAbsorption(double z);
  void computeElectronAbsorption(double z);
};

}  // namespace gnuprop

#endif  // GAMMAPROP_H