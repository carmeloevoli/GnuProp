#ifndef SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H
#define SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H

#include <memory>

namespace beniamino {

// class CosmoNeutrinos {
//  public:
//   CosmoNeutrinos(const solutions::Beniamino& b,
//                  const std::shared_ptr<photonfields::PhotonField>& ebl);
//   virtual ~CosmoNeutrinos() = default;

//   double computeNeutrinoFlux(double Enu, double zMax, size_t N = 9) const;
//   double getProtonFlux(double Ep, double z) const;

//  protected:
//   double I_deps(double EnuObs, double Ep, double z, size_t N = 10) const;
//   double I_dEp(double EnuObs, double z, size_t N = 10) const;

//  protected:
//   std::shared_ptr<photonfields::PhotonField> m_ebl;
//   std::shared_ptr<cosmo::Cosmology> m_cosmology;
//   std::shared_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> m_nuSpec;
//   utils::LookupTable<200, 101> m_Jp;
// };

}  // namespace beniamino

#endif  // SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H