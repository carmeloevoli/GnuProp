#include "cosmoneutrinos.h"

#include "simprop/utils/numeric.h"

namespace beniamino {

CosmoNeutrinos::CosmoNeutrinos(const Beniamino& b) {
  {
    auto f = [&b](double logEp, double z) -> double {
      auto value = b.computeFlux(std::exp(logEp), z, 1e-4);
      assert(value >= 0.);
      return std::log(std::max(value, 1e-30));
    };
    const auto zMax = b.getMaxRedshift();
    m_Jp.cacheTable(f, {log(1e15 * SI::eV), log(1e22 * SI::eV)}, {0., zMax});
  }
}

double CosmoNeutrinos::I_deps(double Enu, double Ep, double z, size_t N) const {
  if (Enu >= Ep) return 0;
  const auto epsMin = m_cmb.getMinPhotonEnergy();
  const auto epsMax = m_cmb.getMaxPhotonEnergy();
  auto integrand = [this, Ep, Enu, z](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    auto value = (eta < 1.) ? 0. : m_cmb.density(epsilon, z) * m_nuSpec.Phi(eta, Enu / Ep);
    return epsilon * value;
  };
  auto a = std::log(epsMin);
  auto b = std::log(epsMax);

  return simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
}

double CosmoNeutrinos::nuEmissivity(double Enu, double z, size_t N) const {
  auto integrand = [this, Enu, z](double lnEp) {
    const auto Ep = std::exp(lnEp);
    auto lgJp = m_Jp.get(lnEp, z);
    if (lgJp < -80) return 0.;
    auto value = std::exp(m_Jp.get(lnEp, z));
    value *= I_deps(Enu, Ep, z);
    return value;
  };
  auto a = std::log(Enu);
  auto b = std::log(std::min(1e5 * Enu, 1e21 * SI::eV));

  auto result = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
  return 4. * M_PI / SI::cLight * result;
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