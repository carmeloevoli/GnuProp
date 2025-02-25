#ifndef GNUPROP_CACHED_GAMMAABSORPTIONCACHED_H
#define GNUPROP_CACHED_GAMMAABSORPTIONCACHED_H

#include "cached/cached.h"
#include "interactions/PhotoPair.h"
#include "simprop.h"

namespace cache {

class GammaAbsorptionCached {
 private:
  const double m_units = 1. / SI::Gyr;
  std::vector<double> m_energyAxis;
  std::vector<double> m_redshiftAxis;
  std::shared_ptr<simprop::photonfields::PhotonField> m_phField;
  std::unique_ptr<Interactions::PhotoPair> m_photoPair =
      std::make_unique<Interactions::PhotoPair>();

 public:
  GammaAbsorptionCached() {};
  virtual ~GammaAbsorptionCached() = default;

  void buildEnergyAxis(double energyMin, double energyMax, size_t energySize) {
    m_energyAxis = simprop::utils::LogAxis(energyMin, energyMax, energySize);
  }

  void buildRedshiftAxis(double zMin, double zMax, size_t zSize) {
    m_redshiftAxis = simprop::utils::LinAxis(zMin, zMax, zSize);
  }

  void buildPhotonField(std::shared_ptr<simprop::photonfields::PhotonField> phField) {
    m_phField = std::move(phField);
  }

  void run(const std::string& filename) {
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
            const size_t N = 10000;
            auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, N, 1e-3);

            return s * m_photoPair->sigma_com(s) * value;
          };

          const auto sMin = 4. * pow2(SI::electronMassC2);
          const auto sMax = 4. * E_gamma * epsMax;
          const size_t N = 10000;
          auto value = simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, N, 1e-3);
          value *= SI::cLight / 8. / pow2(E_gamma);
          value /= std::pow(1. + z, 3.);

          return std::max(value / m_units, 0.);
        },
        m_redshiftAxis, m_energyAxis);

    cache.computeAndSave();
  }
};

// (std::shared_ptr<simprop::photonfields::PhotonField> phField,
//                          const std::vector<double>& redshifts,
//                          const std::vector<double>& energyAxis, const std::string& filename)
//                          {
//   const auto photoPair = Interactions::PhotoPair();

//   CachedFunction2D cache(
//       filename,
//       [&](double z, double E_gamma) {
//         const auto epsMax = phField->getMaxPhotonEnergy();

//         auto integrandOuter = [&](double s) {
//           const auto epsThr = s / 4. / E_gamma;
//           const auto epsMin = std::max(epsThr, phField->getMinPhotonEnergy());

//           auto integrandInner = [&](double lnEpsilon) {
//             const auto epsilon = std::exp(lnEpsilon);
//             return phField->density(epsilon, z) / epsilon;
//           };

//           const auto a = std::log(epsMin);
//           const auto b = std::log(epsMax);
//           const size_t N = 10000;
//           auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, N, 1e-3);

//           return s * photoPair.sigma_com(s) * value;
//         };

//         const auto sMin = 4. * pow2(SI::electronMassC2);
//         const auto sMax = 4. * E_gamma * epsMax;
//         const size_t N = 10000;
//         auto value = simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, N,
//         1e-3); value *= SI::cLight / 8. / pow2(E_gamma); value /= std::pow(1. + z, 3.);

//         return std::max(value / units, 0.);
//       },
//       redshifts, energyAxis);

//   cache.computeAndSave();
// }

}  // namespace cache

#endif