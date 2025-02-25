#ifndef GNUPROP_CACHED_PROTONLOSSESCACHED_H
#define GNUPROP_CACHED_PROTONLOSSESCACHED_H

#include "cached/cached.h"
#include "simprop.h"

namespace cache {

class ProtonLossesCached {
 private:
  const double m_units = 1. / SI::Gyr;
  std::vector<double> m_energyAxis;
  std::unique_ptr<simprop::losses::ContinuousLosses> m_losses;

 public:
  ProtonLossesCached() {};
  virtual ~ProtonLossesCached() = default;

  void buildEnergyAxis(double energyMin, double energyMax, size_t energySize) {
    m_energyAxis = simprop::utils::LogAxis(energyMin, energyMax, energySize);
  }

  void builLosses(std::unique_ptr<simprop::losses::ContinuousLosses> losses) {
    m_losses = std::move(losses);
  }

  void run(const std::string& filename) {
    CachedFunction1D cache(
        filename,
        [&](double energy) {
          const auto gammap = energy / SI::protonMassC2;
          return m_losses->beta(simprop::proton, gammap) / m_units;
        },
        m_energyAxis);

    cache.computeAndSave();
  }
};

}  // namespace cache

#endif