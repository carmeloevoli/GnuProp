#ifndef GNUPROP_CACHED_ABSTRACTCACHED_H
#define GNUPROP_CACHED_ABSTRACTCACHED_H

#include "cached/cached.h"

namespace cache {

class AbstractCached {
 protected:
  const size_t LIMIT = 10000;
  const double m_units = 1. / SI::Gyr;
  double m_precision = 1e-4;
  std::vector<double> m_energyAxis;
  std::vector<double> m_xAxis;
  std::vector<double> m_redshiftAxis;
  std::shared_ptr<simprop::photonfields::PhotonField> m_phField;

 public:
  virtual ~AbstractCached() = default;

  void buildEnergyAxis(double energyMin, double energyMax, size_t energySize) {
    m_energyAxis = simprop::utils::LogAxis(energyMin, energyMax, energySize);
  }

  void buildXAxis(double xMin, double xMax, size_t xSize) {
    m_xAxis = simprop::utils::LogAxis(xMin, xMax, xSize);
  }

  void buildRedshiftAxis(double zMin, double zMax, size_t zSize) {
    m_redshiftAxis = simprop::utils::LinAxis(zMin, zMax, zSize);
  }

  void buildPhotonField(std::shared_ptr<simprop::photonfields::PhotonField> phField) {
    m_phField = std::move(phField);
  }

  void setPrecision(double prec) { m_precision = prec; }

  virtual void run(const std::string& filename) = 0;
};

}  // namespace cache

#endif  // GNUPROP_CACHED_ABSTRACTCACHED_H
