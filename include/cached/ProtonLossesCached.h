#ifndef GNUPROP_CACHED_PROTONLOSSESCACHED_H
#define GNUPROP_CACHED_PROTONLOSSESCACHED_H

#include "cached/cached.h"
#include "simprop.h"

namespace cache {

template <class T>
void ProtonLosses(std::shared_ptr<simprop::photonfields::PhotonField> phField,
                  const std::vector<double>& energyAxis, const std::string& filename) {
  const auto losses = T(phField);
  const auto units = 1. / SI::Gyr;

  CachedFunction1D cache(
      filename,
      [&](double energy) {
        const auto gammap = energy / SI::protonMassC2;
        return losses.beta(simprop::proton, gammap) / units;
      },
      energyAxis);

  cache.computeAndSave();
}

}  // namespace cache

#endif