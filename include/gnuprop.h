#ifndef GNUPROP_H
#define GNUPROP_H

#include <memory>

#include "losses.h"
#include "rates.h"
#include "simprop/core/cosmology.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace gnuprop {

class GnuProp {
 public:
  explicit GnuProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology);
  virtual ~GnuProp() = default;

  double getRedshiftMax() const { return m_zMax; }

  void setRedshiftSize(size_t N) { m_zSize = N; }
  void setRedshiftMax(double zMax) { m_zMax = zMax; }
  void enablePhotoPion() { m_doPhotoPion = true; }
  void disablePhotoPion() { m_doPhotoPion = false; }

  // protons
  void evolve(double zObs);
  void dump(const std::string& filename) const;
  void build();

  // neutrinos
  // double nuInteractionRate(double Enu, double Ep, double z, size_t N = 1000) const;

 protected:
  std::unique_ptr<simprop::cosmo::Cosmology> m_cosmology;
  std::unique_ptr<gnuprop::EnergyLosses> m_losses;
  std::unique_ptr<gnuprop::NeutrinoProductionRate> m_nuProductionRate;
  // std::unique_ptr<simprop::photonfields::CMB> m_cmb;
  // std::unique_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> m_nuSpec;

  double m_sourceComovingEmissivity{4.434e45 * SI::erg / SI::Mpc3 / SI::year};
  double m_injSlope{2.6};
  double m_evolutionIndex{0.0};
  double m_expCutoff{1e21 * SI::eV};
  double m_zMax{2.0};
  double m_energyMin{1e17 * SI::eV};
  double m_energyMax{1e23 * SI::eV};
  size_t m_energySize{192};  // TODO improve
  size_t m_zSize{10000};

  bool m_doPhotoPion = false;

  // energy axis
  std::vector<double> m_eAxis;
  // protons
  std::vector<double> m_np;
  // neutrinos
  std::vector<double> m_qnu;
  std::vector<double> m_nnu;

  double Source(double E, double z) const;
  void computeNuEmissivity(double z);

 protected:  // CN
  std::vector<double> knownTerm;
  std::vector<double> diagonal;
  std::vector<double> upperDiagonal;
  std::vector<double> lowerDiagonal;
};

}  // namespace gnuprop

#endif  // GAMMAPROP_H