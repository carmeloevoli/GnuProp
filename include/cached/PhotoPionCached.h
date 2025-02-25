#ifndef GNUPROP_CACHED_PHOTOPIONCACHED_H
#define GNUPROP_CACHED_PHOTOPIONCACHED_H

#include "cached/cached.h"
#include "interactions/PhotoPion.h"
#include "simprop.h"

namespace cache {

template <typename T>
class PhotoPionCached {
 private:
  const double m_units = 1. / SI::Gyr;
  double m_precision = 1e-3;

  std::vector<double> m_energyAxis;
  std::vector<double> m_xAxis;
  std::vector<double> m_redshiftAxis;
  std::shared_ptr<simprop::photonfields::PhotonField> m_phField;
  std::unique_ptr<T> m_secondarySpectrum = std::make_unique<T>();

 public:
  PhotoPionCached() {};
  virtual ~PhotoPionCached() = default;

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

  void run(const std::string& filename) {
    CachedFunction3D cache(
        filename,
        [&](double z, double x, double E_p) {
          if (x > 1.) return 0.;

          const auto m_pi = SI::pionMassC2;
          const auto m_p = SI::protonMassC2;
          const auto epsTh = (pow2(m_pi) + 2. * m_p * m_pi) / 4. / E_p;
          const auto epsMin = std::max(epsTh, m_phField->getMinPhotonEnergy());
          const auto epsMax = m_phField->getMaxPhotonEnergy();

          if (epsMin > epsMax) return 0.;

          auto integrand = [&](double lnEpsilon) {
            const auto epsilon = std::exp(lnEpsilon);
            const auto eta = 4. * epsilon * E_p / pow2(SI::protonMassC2);
            auto value = m_phField->density(epsilon, z) * m_secondarySpectrum->Phi(eta, x);
            return epsilon * value;
          };

          const auto a = std::log(epsMin);
          const auto b = std::log(epsMax);
          const size_t N = 10000;
          auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, m_precision);
          value /= std::pow(1. + z, 3.);

          return std::max(value / m_units, 0.);
        },
        m_redshiftAxis, m_xAxis, m_energyAxis);

    cache.computeAndSave();
  }
};

}  // namespace cache

#endif