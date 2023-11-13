#ifndef BENIAMINO_BENIAMINO_H
#define BENIAMINO_BENIAMINO_H

#include <memory>

#include "losses.h"
#include "simprop/core/cosmology.h"
#include "simprop/utils/logging.h"

namespace beniamino {

struct SourceParams {
  double injSlope;
  double evolutionIndex;
  double expCutoff;
  double zMax;
};

class Uhecr {
 public:
  Uhecr(std::unique_ptr<simprop::cosmo::Cosmology> cosmology);
  Uhecr(std::unique_ptr<simprop::cosmo::Cosmology> cosmology, const SourceParams& params);
  virtual ~Uhecr() = default;

  double generationEnergy(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double dilationFactor(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double relError = 1e-2) const;
  double getMaxRedshift() const { return m_zMax; }

 protected:
  std::unique_ptr<beniamino::LossesTable> m_losses;
  std::unique_ptr<simprop::cosmo::Cosmology> m_cosmology;

  const double m_sourceEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  const double m_maxEnergy{1e24 * SI::eV};
  const double m_minEnergy{1e17 * SI::eV};

  double m_zMax{3.};
  double m_injSlope{2.6};
  double m_evolutionIndex{0.};
  double m_expCutoff{-1.};
};

}  // namespace beniamino

#endif