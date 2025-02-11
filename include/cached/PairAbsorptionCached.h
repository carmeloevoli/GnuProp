#ifndef GNUPROP_CACHED_PAIRABSORPTIONCACHED_H
#define GNUPROP_CACHED_PAIRABSORPTIONCACHED_H

#include "cached/cached.h"
#include "interactions/InverseCompton.h"
#include "simprop.h"

namespace cache {

void PairAbsorptionRate(std::shared_ptr<simprop::photonfields::PhotonField> phField,
                        const std::vector<double>& redshifts, const std::vector<double>& energyAxis,
                        const std::string& filename) {
  const auto units = 1. / SI::Gyr;
  const auto ic = Interactions::InverseCompton();

  CachedFunction2D cache(
      filename,
      [&](double z, double E_pair) {
        const auto me_2 = pow2(SI::electronMassC2);
        const auto beta_e = std::sqrt(1. - me_2 / pow2(E_pair));
        const auto epsMax = phField->getMaxPhotonEnergy();

        auto integrandOuter = [&](double s) {
          const auto epsThr = (s - me_2) / 2. / E_pair / (1. + beta_e);
          const auto epsMin = std::max(epsThr, phField->getMinPhotonEnergy());

          auto integrandInner = [&](double lnEpsilon) {
            const auto epsilon = std::exp(lnEpsilon);
            return phField->density(epsilon, z) / epsilon;
          };

          const auto a = std::log(epsMin);
          const auto b = std::log(epsMax);
          const size_t N = 10000;
          auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, N, 1e-3);

          return (s - me_2) * ic.sigma_com(s) * value;
        };

        const auto sMin = me_2;
        const auto sMax = me_2 + 2. * E_pair * epsMax * (1. + beta_e);

        const size_t N = 10000;
        auto value = simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, N, 1e-3);
        value *= SI::cLight / 8. / beta_e / pow2(E_pair);
        value /= std::pow(1. + z, 3.);

        return std::max(value / units, 0.);
      },
      redshifts, energyAxis);

  cache.computeAndSave();
}

}  // namespace cache

#endif