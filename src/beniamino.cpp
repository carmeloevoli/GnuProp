#include "beniamino.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "cosmology.h"
#include "numeric.h"

namespace beniamino {

#define INTSTEPS 1000
#define VERYSMALLENERGY (1e15 * SI::eV)
#define VERYLARGEENERGY (1e30 * SI::eV)
#define VERYLARGEJACOBIAN (1e5)

Beniamino::Beniamino(const SourceParams &params) {
  m_injSlope = params.injSlope;
  m_evolutionIndex = params.evolutionIndex;
  m_expCutoff = params.expCutoff;
  m_zMax = params.zMax;

  m_losses = std::make_shared<beniamino::LossesTable<double>>("data/SimProp_proton_losses.txt");
  assert(m_losses->loadTable());
}

double Beniamino::generationEnergy(double E, double zInit, double zFinal, double relError) const {
  assert(zFinal >= zInit);
  assert(E > 0);
  auto dEdz = [this](double z, double E_g) {
    if (E_g < VERYLARGEENERGY) {
      auto E = E_g * (1. + z);
      return E_g * (1. / (1. + z) + dtdz(z) * pow3(1. + z) * m_losses->beta(E));
    } else
      return 0.;
  };
  auto value = utils::odeiv<double>(dEdz, E, zInit, zFinal, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEENERGY);
}

double Beniamino::dilationFactor(double E, double zInit, double zFinal, double relError) const {
  assert(zFinal >= zInit);
  assert(E > 0);
  auto dydz = [this, E, zInit](double z, double y) {
    if (y < VERYLARGEJACOBIAN) {
      auto E_g = generationEnergy(E, zInit, z, 1e-5);
      auto E_prime = std::min(E_g * (1. + z), VERYLARGEENERGY);
      auto dbdE = std::max(m_losses->dbdE(E_prime), 0.);
      return y * (1. / (1. + z) + dtdz(z) * pow3(1. + z) * dbdE);
    } else
      return 0.;
  };
  auto value = utils::odeiv<double>(dydz, 1., zInit, zFinal, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEJACOBIAN);
}

double Beniamino::computeFlux(double E, double zObs, double relError) const {
  if (zObs >= m_zMax) return 0;
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
                                             0.001);  // TODO check this
  return factor * I;
}

}  // namespace beniamino