#include "losses.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "simprop/core/units.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/utils/numeric.h"
#include "utils.h"

namespace beniamino {

EnergyLosses::EnergyLosses() {
  auto cmb = std::make_shared<simprop::photonfields::CMB>();
  m_pair = std::make_shared<simprop::losses::PairProductionLosses>(cmb);
  m_photopi = std::make_shared<simprop::losses::PhotoPionContinuousLosses>(cmb);
}

double EnergyLosses::beta(double E, double z) const {
  return m_pair->beta(simprop::proton, E / SI::protonMassC2, z) +
         m_photopi->beta(simprop::proton, E / SI::protonMassC2, z);
}

// double EnergyLosses::dbdE(double E) const {
//   auto dbetadlnE = simprop::utils::deriv<double>(
//       [this](double lnx) {
//         auto x = std::exp(lnx);
//         return beta(x);
//       },
//       std::log(E), 0.05);
//   return beta(E) + dbetadlnE;
// }

}  // namespace beniamino