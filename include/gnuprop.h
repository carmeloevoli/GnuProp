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

  // secondary models
  void addNeutrinoSource(std::unique_ptr<gnuprop::ProductionRate> q_nu) {
    m_photoPionNus.push_back(std::move(q_nu));
  }
  // secondary models
  void addPhotonSource(std::unique_ptr<gnuprop::ProductionRate> q_gamma) {
    m_photoPionGammas.push_back(std::move(q_gamma));
  }

  // gnuprop
  void build();
  void evolve(double zObs);
  void dump(const std::string& filename) const;

 protected:  // main parameters
  double m_sourceComovingEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  double m_injSlope{2.5};
  double m_evolutionIndex{0.0};
  double m_heCutoff{1e23 * SI::eV};
  double m_leCutoff{1e17 * SI::eV};
  double m_zMax{5.0};
  double m_energyMin{1e15 * SI::eV};
  double m_energyMax{1e22 * SI::eV};
  size_t m_energySize{7 * 64};
  size_t m_zSize{10000};

 protected:
  // cosmology
  std::unique_ptr<simprop::cosmo::Cosmology> m_cosmology;

  // energy axis
  std::vector<double> m_eAxis;

  // protons
  std::vector<std::unique_ptr<gnuprop::ProtonLossRate>> m_losses;
  std::vector<double> m_qp;
  std::vector<double> m_betap;
  std::vector<double> m_np;

  // neutrinos
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_photoPionNus;
  std::vector<double> m_qnu;
  std::vector<double> m_nnu;

  // photons
  std::vector<std::unique_ptr<gnuprop::ProductionRate>> m_photoPionGammas;
  std::unique_ptr<gnuprop::AbsorptionRate> m_absGammas;
  std::vector<double> m_qgamma;
  std::vector<double> m_ngamma;

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
};

}  // namespace gnuprop

#endif  // GAMMAPROP_H