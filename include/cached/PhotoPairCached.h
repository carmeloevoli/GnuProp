#ifndef GNUPROP_CACHED_PHOTOPAIRCACHED_H
#define GNUPROP_CACHED_PHOTOPAIRCACHED_H

#include "cached/cached.h"
#include "interactions/PhotoPair.h"
#include "simprop.h"

namespace cache {

void PhotoPairRate(std::shared_ptr<simprop::photonfields::PhotonField> phField,
                   const std::vector<double>& redshifts, const std::vector<double>& xAxis,
                   const std::vector<double>& eGammaAxis, const std::string& filename) {
  const auto units = 1. / SI::Gyr;
  const auto photoPair = Interactions::PhotoPair();

  CachedFunction3D cache(
      filename,
      [&](double z, double x, double E_gamma) {
        if (x > 1.) return 0.;

        const auto me_2 = pow2(SI::electronMassC2);
        const auto epsTh = me_2 / E_gamma;
        const auto epsMin = std::max(epsTh, phField->getMinPhotonEnergy());
        const auto epsMax = phField->getMaxPhotonEnergy();

        if (epsMin > epsMax) return 0.;

        auto integrand = [&](double lnEpsilon) {
          const auto epsilon = std::exp(lnEpsilon);
          auto value =
              phField->density(epsilon, z) * photoPair.dsigma_dE(E_gamma, epsilon, x * E_gamma);

          return epsilon * value;
        };

        const auto a = std::log(epsMin);
        const auto b = std::log(epsMax);
        const size_t N = 10000;
        auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
        value /= std::pow(1. + z, 3.);

        return std::max(value / units, 0.);
      },
      redshifts, xAxis, eGammaAxis);

  cache.computeAndSave();
}

}  // namespace cache

#endif
