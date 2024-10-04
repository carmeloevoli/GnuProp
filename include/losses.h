#ifndef BENIAMINO_LOSSES_H
#define BENIAMINO_LOSSES_H

#include <string>
#include <vector>

#include "simprop/energyLosses/PairProductionLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

namespace beniamino {

class EnergyLosses {
 public:
  EnergyLosses();
  double beta(double E, double z = 0) const;
  // double dbdE(double E) const;

 private:
  std::shared_ptr<simprop::losses::PairProductionLosses> m_pair;
  std::shared_ptr<simprop::losses::PhotoPionContinuousLosses> m_photopi;
};

}  // namespace beniamino

#endif