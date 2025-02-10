#ifndef GNUPROP_CACHED_PHOTOPIONCACHED_H
#define GNUPROP_CACHED_PHOTOPIONCACHED_H

#include "cached/cached.h"
#include "interactions/PhotoPion.h"
#include "simprop.h"

namespace cache {

template <class T>
void PhotoPionRate(std::shared_ptr<simprop::photonfields::PhotonField> phField,
                   const std::vector<double>& redshifts, const std::vector<double>& xAxis,
                   const std::vector<double>& eProtonAxis, const std::string& filename) {
  const auto units = 1. / SI::Gyr;
  const auto secSpec = T();

  CachedFunction3D cache(
      filename,
      [&](double z, double x, double E_p) {
        if (x > 1.) return 0.;

        const auto m_pi = SI::pionMassC2;
        const auto m_p = SI::protonMassC2;
        const auto epsTh = (pow2(m_pi) + 2. * m_p * m_pi) / 4. / E_p;
        const auto epsMin = std::max(epsTh, phField->getMinPhotonEnergy());
        const auto epsMax = phField->getMaxPhotonEnergy();

        if (epsMin > epsMax) return 0.;

        auto integrand = [&](double lnEpsilon) {
          const auto epsilon = std::exp(lnEpsilon);
          const auto eta = 4. * epsilon * E_p / pow2(SI::protonMassC2);
          auto value = phField->density(epsilon, z) * secSpec.Phi(eta, x);
          return epsilon * value;
        };

        const auto a = std::log(epsMin);
        const auto b = std::log(epsMax);
        const size_t N = 10000;
        auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
        value /= std::pow(1. + z, 3.);

        return std::max(value / units, 0.);
      },
      redshifts, xAxis, eProtonAxis);

  cache.computeAndSave();
}

}  // namespace cache

#endif