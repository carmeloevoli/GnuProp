#ifndef GNUPROP_CACHED_GAMMAABSORPTIONCACHED_H
#define GNUPROP_CACHED_GAMMAABSORPTIONCACHED_H

#include "cached/cached.h"
#include "interactions/GammaPairProduction.h"
#include "simprop.h"

namespace cache {

void GammaAbsorptionRate(std::shared_ptr<simprop::photonfields::PhotonField> phField,
                         const std::vector<double>& redshifts,
                         const std::vector<double>& energyAxis, const std::string& filename) {
  const auto units = 1. / SI::Gyr;
  const auto sigma_pp = Interactions::GammaPairProduction();

  CachedFunction2D cache(
      filename,
      [&](double z, double E_gamma) {
        const auto epsMax = phField->getMaxPhotonEnergy();

        auto integrandOuter = [&](double s) {
          const auto epsThr = s / 4. / E_gamma;
          const auto epsMin = std::max(epsThr, phField->getMinPhotonEnergy());

          auto integrandInner = [&](double lnEpsilon) {
            const auto epsilon = std::exp(lnEpsilon);
            return phField->density(epsilon, z) / epsilon;
          };

          const auto a = std::log(epsMin);
          const auto b = std::log(epsMax);
          const size_t N = 10000;
          auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, N, 1e-3);

          return s * sigma_pp.sigmaInCoMFrame(s) * value;
        };

        const auto sMin = 4. * pow2(SI::electronMassC2);
        const auto sMax = 4. * E_gamma * epsMax;
        const size_t N = 10000;
        auto value = simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, N, 1e-3);
        value *= SI::cLight / 8. / pow2(E_gamma);
        value /= std::pow(1. + z, 3.);

        return std::max(value / units, 0.);
      },
      redshifts, energyAxis);

  cache.computeAndSave();
}

}  // namespace cache

#endif