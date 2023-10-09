#ifndef BENIAMINO_BENIAMINO_H
#define BENIAMINO_BENIAMINO_H

#include <memory>

#include "losses.h"
#include "simprop/core/cosmology.h"

namespace beniamino {

struct SourceParams {
  double injSlope;
  double evolutionIndex;
  double expCutoff;
  double zMax;
};

class Beniamino {
 public:
  Beniamino();
  Beniamino(const SourceParams &params);
  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double dilationFactor(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double relError = 1e-2) const;

 protected:
  std::shared_ptr<beniamino::LossesTable> m_losses;
  simprop::cosmo::Cosmology m_cosmology;

  const double m_sourceEmissivity{3e45 * SI::erg / SI::Mpc3 / SI::year};
  const double m_maxEnergy{1e24 * SI::eV};
  const double m_minEnergy{1e17 * SI::eV};

  double m_zMax{3.};
  double m_injSlope{2.6};
  double m_evolutionIndex{0.};
  double m_expCutoff{-1.};
};

}  // namespace beniamino

#endif