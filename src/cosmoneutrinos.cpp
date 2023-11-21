#include "cosmoneutrinos.h"

#include "simprop/utils/numeric.h"

namespace beniamino {

CosmoNeutrinos::CosmoNeutrinos(const Uhecr& b) {
  {
    auto f = [&b](double lnEp, double z) -> double {
      auto value = b.computeFlux(std::exp(lnEp), z, 1e-4);
      assert(value >= 0.);
      return std::log(std::max(value, 1e-30));
    };
    const auto zMax = b.getMaxRedshift();
    m_Jp.cacheTable(f, {log(1e15 * SI::eV), log(1e22 * SI::eV)}, {0., zMax});
  }
}

double CosmoNeutrinos::interactionRate(double Enu, double Ep, double z, size_t N) const {
  if (Enu >= Ep) return 0;
  const auto epsTh = (pow2(SI::pionMassC2) + 2. * SI::protonMassC2 * SI::pionMassC2) / 4. / Ep;
  const auto epsMin = std::max(epsTh, m_cmb.getMinPhotonEnergy());
  const auto epsMax = m_cmb.getMaxPhotonEnergy();
  auto integrand = [this, Ep, Enu, z](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    auto value = m_cmb.density(epsilon, z) * m_nuSpec.Phi(eta, Enu / Ep);
    return epsilon * value;
  };

  auto a = std::log(epsMin);
  auto b = std::log(epsMax);
  auto value = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
  return value;
}

double CosmoNeutrinos::nuEmissivity(double Enu, double z, size_t N) const {
  auto integrand = [this, Enu, z](double lnEp) {
    const auto Ep = std::exp(lnEp);
    auto lgJp = m_Jp.get(lnEp, z);
    if (lgJp < -80) return 0.;
    auto value = std::exp(m_Jp.get(lnEp, z));
    value *= interactionRate(Enu, Ep, z);
    return value;
  };

  auto a = std::log(Enu);
  auto b = std::log(std::min(1e5 * Enu, 1e22 * SI::eV));
  auto value = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
  return 4. * M_PI / SI::cLight * value;
}

double CosmoNeutrinos::computeNeutrinoFlux(double EnuObs, double zMax, size_t N) const {
  return 0.;
  // const auto factor = 1.;
  // for (size_t i = 1; i < 20; i++) {
  //   auto Enu = 1e16 * SI::eV;
  //   std::cout << i << " " << I_dEp(Enu, 0, i) << "\n";
  // }
  // exit(1);
  // auto integrand = [this, EnuObs](double z) {
  //   return m_cosmology.dtdz(z) / pow2(1. + z) * I_dEp(EnuObs * (1. + z), z);
  // };
  // auto I = simprop::utils::RombergIntegration<double>(integrand, 0., zMax, N, 1e-2);
  // return factor * I;
}

}  // namespace beniamino