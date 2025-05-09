#ifndef GNUPROP_CACHED_GAMMAABSORPTIONCACHED_H
#define GNUPROP_CACHED_GAMMAABSORPTIONCACHED_H

#include "cached/AbstractCached.h"
#include "cached/cached.h"
#include "interactions/PhotoPair.h"
#include "simprop.h"

namespace cache {

class GammaAbsorptionCached : public AbstractCached {
 private:
  std::unique_ptr<Interactions::PhotoPair> m_photoPair =
      std::make_unique<Interactions::PhotoPair>();

 public:
  virtual ~GammaAbsorptionCached() = default;

  void run(const std::string& filename, double precision = 1e-3) override {
    CachedFunction2D cache(
        filename,
        [&](double z, double E_gamma) {
          const auto epsMax = m_phField->getMaxPhotonEnergy();

          auto integrandOuter = [&](double s) {
            const auto epsThr = s / 4. / E_gamma;
            const auto epsMin = std::max(epsThr, m_phField->getMinPhotonEnergy());

            auto integrandInner = [&](double lnEpsilon) {
              const auto epsilon = std::exp(lnEpsilon);
              return m_phField->density(epsilon, z) / epsilon;
            };

            const auto a = std::log(epsMin);
            const auto b = std::log(epsMax);
            auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, LIMIT, 1e-3);

            return s * m_photoPair->sigma_com(s) * value;
          };

          const auto sMin = 4. * pow2(SI::electronMassC2);
          const auto sMax = 4. * E_gamma * epsMax;
          auto value =
              simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, LIMIT, precision);
          value *= SI::cLight / 8. / pow2(E_gamma);
          value /= std::pow(1. + z, 3.);

          return std::max(value / m_units, 0.);
        },
        m_redshiftAxis, m_eAxis);

    cache.computeAndSave();
  }
};

}  // namespace cache

#endif