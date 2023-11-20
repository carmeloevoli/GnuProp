#ifndef BENIAMINO_COSMONUS_H
#define BENIAMINO_COSMONUS_H

#include <memory>

#include "KelnerAharonian2008.h"
#include "simprop/core/cosmology.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/utils/lookupContainers.h"
#include "uhecr.h"

namespace beniamino {

class CosmoNeutrinos {
 public:
  CosmoNeutrinos(const Uhecr& b);
  virtual ~CosmoNeutrinos() = default;
  double computeNeutrinoFlux(double Enu, double zMax, size_t N = 10) const;
  double nuEmissivity(double EnuObs, double z, size_t N = 9) const;

 protected:
  simprop::cosmo::Cosmology m_cosmology;
  simprop::photonfields::CMB m_cmb;
  KelnerAharonian2008::NeutrinoProductionSpectrum m_nuSpec;
  simprop::utils::LookupTable<200, 101> m_Jp;

 protected:
  double I_deps(double EnuObs, double Ep, double z, size_t N = 4) const;
};

}  // namespace beniamino

#endif  // SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H