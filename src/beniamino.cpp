#include "beniamino.h"

#include <cassert>

#include "simprop/utils/numeric.h"

namespace beniamino {

#define VERYLARGEENERGY (1e30 * SI::eV)
#define VERYLARGEJACOBIAN (1e6)

Beniamino::Beniamino(std::unique_ptr<simprop::cosmo::Cosmology> cosmology)
    : m_cosmology(std::move(cosmology)) {
  LOGI << "h : " << m_cosmology->h;
  LOGI << "Omega_M : " << m_cosmology->OmegaM;
  LOGI << "Omega_L : " << m_cosmology->OmegaL;
  m_losses = std::make_unique<beniamino::LossesTable>();
  assert(m_losses->loadTable("data/SimProp_proton_losses.txt"));
}

// Beniamino::Beniamino(const SourceParams &params) : Beniamino() {
//   m_injSlope = params.injSlope;
//   m_evolutionIndex = params.evolutionIndex;
//   m_expCutoff = params.expCutoff;
//   m_zMax = params.zMax;
// }

double Beniamino::generationEnergy(double E, double zInit, double zFinal, double relError) const {
  assert(zFinal >= zInit);
  assert(E > 0);
  auto dEdz = [this](double z, double E_g) {
    if (E_g < VERYLARGEENERGY) {
      auto E = E_g * (1. + z);
      return E_g * (1. / (1. + z) + m_cosmology->dtdz(z) * pow3(1. + z) * m_losses->beta(E));
    } else
      return 0.;
  };
  auto value = simprop::utils::odeiv<double>(dEdz, E, zInit, zFinal, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEENERGY);
}

double Beniamino::dilationFactor(double E, double zInit, double zFinal, double relError) const {
  assert(zFinal >= zInit);
  assert(E > 0);
  auto dydz = [this, E, zInit, relError](double z, double y) {
    if (y < VERYLARGEJACOBIAN) {
      auto E_g = generationEnergy(E, zInit, z, 0.1 * relError);
      auto E_prime = std::min(E_g * (1. + z), VERYLARGEENERGY);
      auto dbdE = std::max(m_losses->dbdE(E_prime), 0.);
      return y * (1. / (1. + z) + m_cosmology->dtdz(z) * pow3(1. + z) * dbdE);
    } else
      return 0.;
  };
  auto value = simprop::utils::odeiv<double>(dydz, 1., zInit, zFinal, relError);
  assert(value > 0.);

  return std::min(value, VERYLARGEJACOBIAN);
}

double Beniamino::computeFlux(double E, double zObs, double relError) const {
  if (zObs >= m_zMax) return 0;
  const auto K = (m_injSlope - 2.) / pow2(m_minEnergy);
  const auto L_0 = m_sourceEmissivity;
  const auto factor = SI::cLight / 4. / M_PI * K * L_0;
  auto integrand = [this, E, zObs](double z) {
    const auto E_g = generationEnergy(E, zObs, z, 1e-5);
    if (E_g > m_maxEnergy) return 0.;
    auto dEgdE = dilationFactor(E, zObs, z, 1e-5);
    auto inj = std::pow(E_g / m_minEnergy, -m_injSlope);
    if (m_expCutoff > 0.) inj *= std::exp(-E_g / m_expCutoff);
    auto sourceEvolution = std::pow(1. + z, m_evolutionIndex);
    return m_cosmology->dtdz(z) * sourceEvolution * inj * dEgdE;
  };
  auto I = simprop::utils::RombergIntegration<double>(integrand, zObs, m_zMax, 15,
                                                      1e-4);  // TODO check this
  return factor * I;
}

}  // namespace beniamino