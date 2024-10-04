#ifndef GAMMAPROP_H
#define GAMMAPROP_H

#include <memory>
#include <vector>

#include "KelnerAharonian2008.h"
#include "losses.h"
#include "simprop/core/cosmology.h"
#include "simprop/photonfields/CmbPhotonField.h"
#include "simprop/utils/logging.h"

namespace gammaprop {

class GammaProp {
 public:
  // Constructor with cosmology only
  explicit GammaProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology);

  // Default virtual destructor
  virtual ~GammaProp() = default;

  // Getter for maximum redshift
  double getMaxRedshift() const { return m_zMax; }

  void setSizeRedshift(size_t N) { m_zSize = N; }

  void evolve(double zObs);
  double pComovingSource(double E, double z) const;
  std::vector<double> computeNuEmissivity(double z) const;
  double nuInteractionRate(double Enu, double Ep, double z, size_t N = 5) const;

  void dump(std::string filename) const;

 protected:
  std::unique_ptr<simprop::cosmo::Cosmology> m_cosmology;
  std::unique_ptr<beniamino::EnergyLosses> m_losses;
  std::unique_ptr<simprop::photonfields::CMB> m_cmb;
  std::unique_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> m_nuSpec;

  double m_sourceEmissivity{4.434e45 * SI::erg / SI::Mpc3 / SI::year};
  double m_injSlope{2.6};
  double m_evolutionIndex{0.0};
  double m_expCutoff{1e22 * SI::eV};
  double m_zMax{5.0};
  double m_energyMin{1e17 * SI::eV};
  double m_energyMax{1e23 * SI::eV};
  size_t m_energySize{400};
  size_t m_zSize{100000};

  std::vector<double> m_eAxis;
  std::vector<double> m_np;
  std::vector<double> m_nnu;
  std::vector<double> m_Qnu;

  void build();
};

}  // namespace gammaprop

#endif  // GAMMAPROP_H