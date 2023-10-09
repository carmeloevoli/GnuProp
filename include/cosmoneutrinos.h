#ifndef BENIAMINO_COSMONUS_H
#define BENIAMINO_COSMONUS_H

#include <memory>

#include "beniamino.h"
#include "simprop/core/cosmology.h"

namespace beniamino {

class CosmoNeutrinos {
 public:
  CosmoNeutrinos(const Beniamino& b);
  virtual ~CosmoNeutrinos() = default;

  double computeNeutrinoFlux(double Enu, double zMax, size_t N = 10) const;

  // double getProtonFlux(double Ep, double z) const;

 protected:
  std::shared_ptr<simprop::cosmo::Planck2018> m_cosmology;

  double I_deps(double EnuObs, double Ep, double z, size_t N = 10) const;
  double I_dEp(double EnuObs, double z, size_t N = 10) const;

  //  protected:
  //   std::shared_ptr<photonfields::PhotonField> m_ebl;
  //   std::shared_ptr<cosmo::Cosmology> m_cosmology;
  //   std::shared_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> m_nuSpec;
  //   utils::LookupTable<200, 101> m_Jp;
};

}  // namespace beniamino

#endif  // SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H