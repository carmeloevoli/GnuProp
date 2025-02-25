#ifndef GNUPROP_CACHED_PAIRABSORPTIONCACHED_H
#define GNUPROP_CACHED_PAIRABSORPTIONCACHED_H

#include "cached/cached.h"
#include "interactions/InverseCompton.h"
#include "simprop.h"

namespace cache {

class PairAbsorptionCached {
 private:
  const double m_units = 1. / SI::Gyr;
  std::vector<double> m_energyAxis;
  std::vector<double> m_redshiftAxis;
  std::shared_ptr<simprop::photonfields::PhotonField> m_phField;
  std::unique_ptr<Interactions::InverseCompton> m_ic =
      std::make_unique<Interactions::InverseCompton>();

 public:
  PairAbsorptionCached() {};
  virtual ~PairAbsorptionCached() = default;

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
        [&](double z, double E_pair) {
          const auto me_2 = pow2(SI::electronMassC2);
          const auto beta_e = std::sqrt(1. - me_2 / pow2(E_pair));
          const auto epsMax = m_phField->getMaxPhotonEnergy();

          auto integrandOuter = [&](double s) {
            const auto epsThr = (s - me_2) / 2. / E_pair / (1. + beta_e);
            const auto epsMin = std::max(epsThr, m_phField->getMinPhotonEnergy());

            auto integrandInner = [&](double lnEpsilon) {
              const auto epsilon = std::exp(lnEpsilon);
              return m_phField->density(epsilon, z) / epsilon;
            };

            const auto a = std::log(epsMin);
            const auto b = std::log(epsMax);
            const size_t N = 10000;
            auto value = simprop::utils::QAGIntegration<double>(integrandInner, a, b, N, 1e-3);

            return (s - me_2) * m_ic->sigma_com(s) * value;
          };

          const auto sMin = me_2;
          const auto sMax = me_2 + 2. * E_pair * epsMax * (1. + beta_e);

          const size_t N = 10000;
          auto value = simprop::utils::QAGIntegration<double>(integrandOuter, sMin, sMax, N, 1e-3);
          value *= SI::cLight / 8. / beta_e / pow2(E_pair);
          value /= std::pow(1. + z, 3.);

          return std::max(value / m_units, 0.);
        },
        m_redshiftAxis, m_energyAxis);

    cache.computeAndSave();
  }
};

}  // namespace cache

#endif