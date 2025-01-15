#include "losses.h"

#include "simprop/core/units.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace gnuprop {

EnergyLosses::EnergyLosses()
    : m_pair(std::make_shared<simprop::photonfields::CMB>()),
      m_photopi(std::make_shared<simprop::photonfields::CMB>()) {}

double EnergyLosses::beta(double E, double z) const {
  const double gammap = E / SI::protonMassC2;
  return m_pair.beta(simprop::proton, gammap, z) + m_photopi.beta(simprop::proton, gammap, z);
}

}  // namespace gnuprop

// double EnergyLosses::dbdE(double E) const {
//   auto dbetadlnE = simprop::utils::deriv<double>(
//       [this](double lnx) {
//         auto x = std::exp(lnx);
//         return beta(x);
//       },
//       std::log(E), 0.05);
//   return beta(E) + dbetadlnE;
// }