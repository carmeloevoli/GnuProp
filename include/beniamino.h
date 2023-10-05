#ifndef BENIAMINO_BENIAMINO_H
#define BENIAMINO_BENIAMINO_H

// #include "simprop/core/cosmology.h"
// #include "simprop/energyLosses/ContinuousLosses.h"
// #include "simprop/utils/lookupContainers.h"

#include <memory>

#include "losses.h"
#include "units.h"

namespace beniamino {

struct SourceParams {
  double injSlope;
  double evolutionIndex;
  double expCutoff;
  double zMax;
};

class Beniamino {
 public:
  Beniamino(const SourceParams &params);
  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double dilationFactor(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double relError = 1e-2) const;

 protected:
  std::shared_ptr<beniamino::LossesTable<double>> m_losses;

  const double m_sourceEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  const double m_maxEnergy{1e24 * SI::eV};
  const double m_minEnergy{1e17 * SI::eV};

  double m_zMax{3.};
  double m_injSlope{2.6};
  double m_evolutionIndex{0.};
  double m_expCutoff{-1.};
};

// class Beniamino {
// public:
//   Beniamino(
//       const SourceParams &params,
//       const std::shared_ptr<cosmo::Cosmology> &cosmology,
//       const std::vector<std::shared_ptr<losses::ContinuousLosses>> &losses);
//   Beniamino &doCaching();

//   virtual ~Beniamino() = default;

//   double dilationFactor(double E, double zInit, double zFinal,
//                         double relError = 1e-3) const;
//
//   // double computeFluxUnm(double E, double zMax, double relError = 1e-3)
//   const;

// public:
//   double getMaxRedshift() const { return m_zMax; }
//   const std::shared_ptr<cosmo::Cosmology> &getCosmology() const {
//     return m_cosmology;
//   }

// protected:
//   std::shared_ptr<cosmo::Cosmology> m_cosmology;
//   std::vector<std::shared_ptr<losses::ContinuousLosses>> m_losses;
//   utils::LookupArray<10000> m_lossesLookup;

//   double m_zMax{3.};
//   double m_injSlope{2.6};
//   double m_evolutionIndex{0.};
//   double m_expCutoff{-1.};
//   bool m_doCaching{false};
// };

}  // namespace beniamino

#endif