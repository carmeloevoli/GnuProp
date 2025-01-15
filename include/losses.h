#ifndef GNUPROP_LOSSES_H
#define GNUPROP_LOSSES_H

#include "simprop/energyLosses/PairProductionLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"

namespace gnuprop {

class EnergyLosses {
 public:
  EnergyLosses();

  // Compute beta
  double beta(double E, double z = 0) const;

  // PhotoPion control
  void enablePhotoPion() { doPhotoPion = true; }
  void disablePhotoPion() { doPhotoPion = false; }
  bool isPhotoPionEnabled() const { return doPhotoPion; }

 private:
  bool doPhotoPion = false;  // PhotoPion flag

  simprop::losses::PairProductionLosses m_pair;
  simprop::losses::PhotoPionContinuousLosses m_photopi;
};

}  // namespace gnuprop

#endif  // GNUPROP_LOSSES_H
