#ifndef GNUPROP_CACHED_ABSTRACTCACHED_H
#define GNUPROP_CACHED_ABSTRACTCACHED_H

#include "cached/cached.h"

namespace cache {

class AbstractCached {
 protected:
  const size_t LIMIT = 10000;
  const double m_units = 1. / SI::Gyr;
  std::vector<double> m_primaryAxis;
  std::vector<double> m_secondaryAxis;
  std::vector<double> m_redshiftAxis;
  std::shared_ptr<simprop::photonfields::PhotonField> m_phField;

 public:
  virtual ~AbstractCached() = default;

  void buildPrimaryEnergyAxis(double energyMin, double energyMax, size_t energySize) {
    m_primaryAxis = simprop::utils::LogAxis(energyMin, energyMax, energySize);
  }

  void buildSecondaryEnergyAxis(double energyMin, double energyMax, size_t energySize) {
    m_secondaryAxis = simprop::utils::LogAxis(energyMin, energyMax, energySize);
  }

  void buildRedshiftAxis(double zMin, double zMax, size_t zSize) {
    m_redshiftAxis = simprop::utils::LinAxis(zMin, zMax, zSize);
  }

  void buildPhotonField(std::shared_ptr<simprop::photonfields::PhotonField> phField) {
    m_phField = std::move(phField);
  }

  virtual void run(const std::string& filename, double precision = 1e-3) = 0;
};

}  // namespace cache

#endif  // GNUPROP_CACHED_ABSTRACTCACHED_H
