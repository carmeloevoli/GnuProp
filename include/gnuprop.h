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
  void setCutoffEnergy(double Ec) { m_expCutoff = Ec; }
  void setInjectionSlope(double slope) { m_injSlope = slope; }

  // gnuprop
  void build();
  void evolve(double zObs);
  void dump(const std::string& filename) const;

 protected:  // main parameters
  double m_sourceComovingEmissivity{2e45 * SI::erg / SI::Mpc3 / SI::year};
  double m_injSlope{2.5};
  double m_evolutionIndex{0.0};
  double m_expCutoff{1e23 * SI::eV};
  double m_zMax{6.0};
  double m_energyMin{1e17 * SI::eV};
  double m_energyMax{1e22 * SI::eV};
  size_t m_energySize{6 * 64};
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
  std::unique_ptr<gnuprop::PhotoPionRate> m_photoPionNu;
  std::vector<double> m_qnu;
  std::vector<double> m_nnu;

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
};

}  // namespace gnuprop

#endif  // GAMMAPROP_H