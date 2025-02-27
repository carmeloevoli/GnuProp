#ifndef GNUPROP_CACHED_ICCACHED_H
#define GNUPROP_CACHED_ICCACHED_H

#include "cached/cached.h"
#include "interactions/InverseCompton.h"
#include "simprop.h"

namespace cache {

class InverseComptonCached {
 private:
  const double m_units = 1. / SI::Gyr;
  double m_precision = 1e-3;
  bool m_doGammas = true;

  std::vector<double> m_energyAxis;
  std::vector<double> m_xAxis;
  std::vector<double> m_redshiftAxis;
  std::shared_ptr<simprop::photonfields::PhotonField> m_phField;

  std::unique_ptr<Interactions::InverseCompton> m_ic =
      std::make_unique<Interactions::InverseCompton>();

 public:
  InverseComptonCached() {};
  virtual ~InverseComptonCached() = default;

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

  void setDoGammas(bool doGammas) { m_doGammas = doGammas; }

  void run(const std::string& filename) {
    CachedFunction3D cache(
        filename,
        [&](double z, double x, double E_electron) {
          if (x > 1.) return 0.;

          const auto epsTh = 0.;
          const auto epsMin = std::max(epsTh, m_phField->getMinPhotonEnergy());
          const auto epsMax = m_phField->getMaxPhotonEnergy();

          if (epsMin > epsMax) return 0.;

          auto integrand = [&](double lnEpsilon) {
            const auto epsilon = std::exp(lnEpsilon);
            const auto dsigma_dE =
                (m_doGammas) ? m_ic->dsigma_dE(E_electron, epsilon, x * E_electron)
                             : m_ic->dsigma_dE(E_electron, epsilon, (1. - x) * E_electron);
            auto value = m_phField->density(epsilon, z) * dsigma_dE;

            return epsilon * value;
          };

          const auto a = std::log(epsMin);
          const auto b = std::log(epsMax);
          const size_t N = 10000;
          auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, m_precision);
          value *= SI::cLight * E_electron;
          value /= std::pow(1. + z, 3.);

          return std::max(value / m_units, 0.);
        },
        m_redshiftAxis, m_xAxis, m_energyAxis);

    cache.computeAndSave();
  }
};

// void InverseComptonRate(std::shared_ptr<simprop::photonfields::PhotonField> phField,
//                         const std::vector<double>& redshifts, const std::vector<double>& xAxis,
//                         const std::vector<double>& eElectronAxis, const std::string& filename,
//                         bool doGammas = true) {
//   const auto units = 1. / SI::Gyr;
//   const auto ic = Interactions::InverseCompton();

//   CachedFunction3D cache(
//       filename,
//       [&](double z, double x, double E_electron) {
//         if (x > 1.) return 0.;

//         const auto epsTh = 0.;
//         const auto epsMin = std::max(epsTh, phField->getMinPhotonEnergy());
//         const auto epsMax = phField->getMaxPhotonEnergy();

//         if (epsMin > epsMax) return 0.;

//         auto integrand = [&](double lnEpsilon) {
//           const auto epsilon = std::exp(lnEpsilon);
//           const auto dsigma_dE = (doGammas)
//                                      ? ic.dsigma_dE(E_electron, epsilon, x * E_electron)
//                                      : ic.dsigma_dE(E_electron, epsilon, (1. - x) * E_electron);
//           auto value = phField->density(epsilon, z) * dsigma_dE;

//           return epsilon * value;
//         };

//         const auto a = std::log(epsMin);
//         const auto b = std::log(epsMax);
//         const size_t N = 10000;
//         auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
//         value *= SI::cLight * E_electron;
//         value /= std::pow(1. + z, 3.);

//         return std::max(value / units, 0.);
//       },
//       redshifts, xAxis, eElectronAxis);

//   cache.computeAndSave();
// }

}  // namespace cache

#endif
