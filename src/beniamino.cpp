#include "beniamino.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "cosmology.h"
#include "numeric.h"

namespace beniamino {

#define INTSTEPS 1000
#define VERYSMALLENERGY (1e15 * SI::eV)
#define VERYLARGEENERGY (1e25 * SI::eV)
#define VERYLARGEJACOBIAN (1e8)

Beniamino::Beniamino(const SourceParams &params) {
  m_injSlope = params.injSlope;
  m_evolutionIndex = params.evolutionIndex;
  m_expCutoff = params.expCutoff;
  m_zMax = params.zMax;
}

// Beniamino &Beniamino::doCaching() {
//   LOGD << "Start Beniamino chaching...";
//   m_lossesLookup.cacheTable(
//       [this](double lnE) {
//         const auto Gamma = std::exp(lnE) / SI::protonMassC2;
//         return std::accumulate(
//             m_losses.begin(), m_losses.end(), 0.,
//             [Gamma](double sum,
//                     const std::shared_ptr<losses::ContinuousLosses> &l) {
//               return sum + l->beta(proton, Gamma);
//             });
//       },
//       {std::log(VERYSMALLENERGY), std::log(VERYLARGEENERGY)});
//   m_doCaching = true;
//   return *this;
// }

double Beniamino::beta(double E) const {
  if (E < VERYSMALLENERGY || E > VERYLARGEENERGY) return 0;
  return 0;  // TODO implement this
}

double Beniamino::dbdE(double E) const {
  if (E < VERYSMALLENERGY || E > VERYLARGEENERGY) return 0;
  // auto dbetadlnE = utils::deriv<double>(
  //     [this](double lnx) {
  //       auto x = std::exp(lnx);
  //       return beta(x);
  //     },
  //     std::log(E), 1e-2);
  // return beta(E) + dbetadlnE;
  return 0;
}

double Beniamino::generationEnergy(double E, double zInit, double zFinal, double relError) const {
  assert(zFinal >= zInit);
  assert(E > 0);
  auto dEdz = [this](double z, double E_g) {
    auto E = std::min(E_g * (1. + z), VERYLARGEENERGY);
    return E_g * (1. / (1. + z) + dtdz(z) * pow3(1. + z) * beta(E));
  };
  auto value = utils::odeiv<double>(dEdz, E, zInit, zFinal, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEENERGY);
}

double Beniamino::dilationFactor(double E, double zInit, double zFinal, double relError) const {
  assert(zFinal >= zInit);
  assert(E > 0);
  if (generationEnergy(E, zInit, zFinal, 1e-2) > 0.2 * VERYLARGEENERGY) return VERYLARGEJACOBIAN;
  auto dydz = [this, E, zInit](double z, double y) {
    auto E_g = generationEnergy(E, zInit, z, 1e-5);
    auto E_prime = std::min(E_g * (1. + z), VERYLARGEENERGY);
    auto dbdE = std::max(this->dbdE(E_prime), 0.);
    return y * (1. / (1. + z) + dtdz(z) * pow3(1. + z) * dbdE);
  };
  auto value = utils::odeiv<double>(dydz, 1., zInit, zFinal, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEJACOBIAN);
}

double Beniamino::computeFlux(double E, double zObs, double relError) const {
  const auto K = (m_injSlope - 2.) / pow2(m_minEnergy);
  const auto L_0 = m_sourceEmissivity;
  const auto factor = SI::cLight / 4. / M_PI * K * L_0;
  auto integrand = [this, E, zObs](double z) {
    const auto E_g = generationEnergy(E, zObs, z, 1e-4);
    if (E_g > m_maxEnergy) return 0.;
    auto dEgdE = dilationFactor(E, zObs, z, 1e-4);
    auto inj = std::pow(E_g / m_minEnergy, -m_injSlope);
    if (m_expCutoff > 0.) inj *= std::exp(-E_g / m_expCutoff);
    auto sourceEvolution = std::pow(1. + z, m_evolutionIndex);
    return dtdz(z) * sourceEvolution * inj * dEgdE;
  };
  auto I = utils::RombergIntegration<double>(integrand, zObs, m_zMax, 8,
                                             0.01);  // TODO check this
  return factor * I;
}

}  // namespace beniamino