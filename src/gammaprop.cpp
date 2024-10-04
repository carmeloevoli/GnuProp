#include "gammaProp.h"

#include <numeric>

namespace gammaprop {

// Constructor with cosmology only
GammaProp::GammaProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology)
    : m_cosmology(std::move(cosmology)) {
  LOGI << "h: " << m_cosmology->h;
  LOGI << "Omega_M: " << m_cosmology->OmegaM;
  LOGI << "Omega_L: " << m_cosmology->OmegaL;
  build();
}

void GammaProp::build() {
  m_eAxis = simprop::utils::LogAxis(m_energyMin, m_energyMax, m_energySize);
  m_np = std::vector<double>(m_energySize, 0.0);
  m_nnu = std::vector<double>(m_energySize, 0.0);
  m_Qnu = std::vector<double>(m_energySize, 0.0);
  m_losses = std::make_unique<beniamino::EnergyLosses>();
  m_cmb = std::make_unique<simprop::photonfields::CMB>();
  m_nuSpec = std::make_unique<KelnerAharonian2008::NeutrinoProductionSpectrum>();
}

void GammaProp::evolve(double zObs) {
  assert(zObs < m_zMax);
  auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
  std::reverse(zAxis.begin(), zAxis.end());
  const auto dz = zAxis[0] - zAxis[1];

  for (const auto& z : zAxis) {
    LOGD << std::setprecision(5) << z << "\t" << std::accumulate(m_np.begin(), m_np.end(), 0.);

    const auto dtdz = std::abs(m_cosmology->dtdz(z));
    const auto H = std::abs(m_cosmology->hubbleRate(z));

    // evolve protons
    std::vector<double> npUp(m_energySize);

    for (size_t i = 0; i < m_energySize - 1.; ++i) {
      const auto E = m_eAxis.at(i);
      const auto Eup = m_eAxis.at(i + 1);
      const auto Q = pow3(1. + z) * pComovingSource(E, z);
      const auto b = E * (m_losses->beta(E, z) + H);
      const auto bUp = 0.;    // Eup * (m_losses->beta(Eup, z) + H);
      const auto dbndE = 0.;  // (bUp * m_np[i + 1] - b * m_np[i]) / (Eup - E);
      const auto value = m_np[i] + dz * (dtdz * Q - 3. / (1. + z) * m_np[i] + dtdz * dbndE);
      assert(value > 0.);
      npUp[i] = value;
    }

    m_np = std::move(npUp);

    // evolve neutrinos
    std::vector<double> nnuUp(m_energySize);

    // m_Qnu = computeNuEmissivity(z);

    for (size_t i = 0; i < m_energySize - 1.; ++i) {
      const auto E = m_eAxis.at(i);
      const auto Eup = m_eAxis.at(i + 1);
      const auto b = E * H;
      const auto bUp = Eup * H;
      const auto dbndE = (bUp * m_nnu[i + 1] - b * m_nnu[i]) / (Eup - E);
      nnuUp[i] = m_nnu[i] + dz * (dtdz * m_Qnu[i] - 3. / (1. + z) * m_nnu[i] + dtdz * dbndE);
    }

    m_nnu = std::move(nnuUp);
  }
}

double GammaProp::pComovingSource(double E, double z) const {
  double value = m_sourceEmissivity * (m_injSlope - 2.0) / (m_energyMin * m_energyMin);
  value *= std::pow(E / m_energyMin, -m_injSlope);
  if (m_expCutoff > 0.0) value *= std::exp(-E / m_expCutoff);
  value *= std::pow(1.0 + z, m_evolutionIndex);
  return std::max(value, 0.0);
}

std::vector<double> GammaProp::computeNuEmissivity(double z) const {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  std::vector<double> Q;
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto Enu = m_eAxis[i];
    double value = 0.;
    for (size_t j = i; j < m_energySize; ++j) {
      value += m_np[j] * nuInteractionRate(Enu, m_eAxis[j], z);
    }
    Q.emplace_back(ln_eRatio * value);
  }

  return Q;
}

double GammaProp::nuInteractionRate(double Enu, double Ep, double z, size_t N) const {
  if (Enu >= Ep) return 0;

  const auto epsTh = (pow2(SI::pionMassC2) + 2. * SI::protonMassC2 * SI::pionMassC2) / 4. / Ep;
  const auto epsMin = std::max(epsTh, m_cmb->getMinPhotonEnergy());
  const auto epsMax = m_cmb->getMaxPhotonEnergy();

  auto integrand = [this, Ep, Enu, z](double lnEpsilon) {
    const auto epsilon = std::exp(lnEpsilon);
    const auto eta = 4. * epsilon * Ep / pow2(SI::protonMassC2);
    auto value = m_cmb->density(epsilon, z) * m_nuSpec->Phi(eta, Enu / Ep);
    return epsilon * value;
  };

  const auto a = std::log(epsMin);
  const auto b = std::log(epsMax);
  // double value = simprop::utils::RombergIntegration<double>(integrand, a, b, N, 1e-3);
  const auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, 1000, 1e-3);
  return value;
}

void GammaProp::dump(std::string filename) const {
  std::ofstream out("output/" + filename);
  const double units = 1. / SI::eV / SI::m2 / SI::sr / SI::sec;
  const double units_emissivity = 1. / SI::eV / SI::m3 / SI::sec;
  out << "# E [eV] - I [eV-1 m-2 s-1 sr-1]\n";
  for (size_t i = 0; i < m_energySize; ++i) {
    out << std::scientific << m_eAxis.at(i) / SI::eV << " ";
    out << SI::cLight / 4. / M_PI * m_np.at(i) / units << " ";
    out << SI::cLight / 4. / M_PI * m_nnu.at(i) / units << " ";
    out << m_Qnu.at(i) / units_emissivity << " ";
    out << "\n";
  }
}

}  // namespace gammaprop
