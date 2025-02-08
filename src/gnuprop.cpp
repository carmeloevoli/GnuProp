#include "gnuprop.h"

#include <cmath>
#include <fstream>
#include <numeric>

#include "simprop.h"
#include "tridiag.h"

namespace gnuprop {

// Constructor with cosmology only
GnuProp::GnuProp(std::unique_ptr<simprop::cosmo::Cosmology> cosmology)
    : m_cosmology(std::move(cosmology)) {
  LOGI << "h: " << m_cosmology->h << " Omega_M: " << m_cosmology->OmegaM
       << " Omega_L: " << m_cosmology->OmegaL;
}

void GnuProp::build() {
  // m_nuProductionRate = std::make_unique<gnuprop::NeutrinoProductionRate>();
  m_losses.push_back(
      std::make_unique<gnuprop::ProtonLossRate>("data/gnuprop_proton_losses_pair.bin"));
  m_losses.push_back(
      std::make_unique<gnuprop::ProtonLossRate>("data/gnuprop_proton_losses_photopion.bin"));
  m_eAxis = simprop::utils::LogAxis(m_energyMin, m_energyMax, m_energySize);

  m_np.assign(m_energySize, 0.0);
  m_betap.assign(m_energySize, 0.0);
  m_qp.assign(m_energySize, 0.0);

  m_nnu.assign(m_energySize, 0.0);
  m_qnu.assign(m_energySize, 0.0);

  knownTerm.assign(m_energySize - 1, 0.0);
  diagonal.assign(m_energySize - 1, 0.0);
  upperDiagonal.assign(m_energySize - 2, 0.0);
  lowerDiagonal.assign(m_energySize - 2, 0.0);
}

void GnuProp::evolve(double zObs) {
  assert(zObs < m_zMax);
  auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
  const auto dz = zAxis[1] - zAxis[0];

  for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {
    auto z = *it;
    // LOGD << std::setprecision(5) << z << "\t" << std::accumulate(m_np.begin(), m_np.end(), 0.);

    const auto dtdz = std::abs(m_cosmology->dtdz(z));
    const auto H = std::abs(m_cosmology->hubbleRate(z));

    // evolve protons
    {
      evolveProtonEmissivity(z);
      evolveProtonLosses(z);
      std::vector<double> npUp(m_energySize - 1, 0.);
      for (size_t i = 0; i < m_energySize - 1; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto dE = Eup - E;
        const auto Q = dtdz * std::pow(1. + z, 3.0) * m_qp[i] - 3. / (1. + z) * m_np[i];
        const auto b = dtdz * E * (m_betap[i] + H);
        const auto bUp = dtdz * Eup * (m_betap[i + 1] + H);
        const auto U_i = 0.5 * dz * bUp / dE;
        const auto C_i = 0.5 * dz * b / dE;
        diagonal[i] = 1. + C_i;
        if (i != m_energySize - 2) upperDiagonal[i] = -U_i;
        knownTerm[i] = U_i * m_np[i + 1] + (1. - C_i) * m_np[i] + dz * Q;
      }

      gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, npUp);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        m_np[i] = std::max(npUp[i], 0.);
      }
    }

    // evolve neutrinos
    {
      evolveNuEmissivity(z);
      std::vector<double> nnuUp(m_energySize, 0.);
      for (size_t i = 0; i < m_energySize - 1.; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto b = E * H;
        const auto bUp = Eup * H;
        const auto dbndE = (bUp * m_nnu[i + 1] - b * m_nnu[i]) / (Eup - E);
        nnuUp[i] = m_nnu[i] + dz * (dtdz * m_qnu[i] - 3. / (1. + z) * m_nnu[i] + dtdz * dbndE);
      }
      m_nnu = std::move(nnuUp);
    }
  }
}

void GnuProp::evolveProtonEmissivity(double z) {
  const auto K = m_sourceComovingEmissivity * (m_injSlope - 2.0) / std::pow(m_energyMin, 2.);
  const auto zfactor = (z <= m_zMax) ? std::pow(1.0 + z, m_evolutionIndex) : 0.;
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E = m_eAxis[i];
    auto value = (E >= m_energyMin) ? std::pow(E / m_energyMin, -m_injSlope) : 0.;
    value *= (m_expCutoff > 0.0) ? std::exp(-E / m_expCutoff) : 1.;
    m_qp[i] = K * zfactor * std::max(value, 0.0);
  }
}

void GnuProp::evolveProtonLosses(double z) {
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto E = m_eAxis[i];
    auto value = 0.;
    for (const auto& loss : m_losses) value += loss->beta(E, z);
    m_betap[i] = value;
  }
}

void GnuProp::evolveNuEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto Enu = m_eAxis[i];
    double value = 0.;
    for (size_t j = i; j < m_energySize; ++j) {
      value += m_np[j] * Enu;  // * m_nuProductionRate->get(Enu, m_eAxis[j], z);
    }
    m_qnu[i] = ln_eRatio * value;
  }
}

void GnuProp::dump(const std::string& filename) const {
  std::ofstream out("output/" + filename);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open output file: output/" + filename);
  }

  const double Ip_units = 1.0 / (SI::eV * SI::m2 * SI::sr * SI::sec);
  const double E2Inu_units = SI::GeV / (SI::cm2 * SI::sr * SI::sec);
  const double q_units = 1.0 / (SI::eV * SI::m3 * SI::sec);
  out << "# E [eV] - I [eV-1 m-2 s-1 sr-1]\n";

  for (size_t i = 0; i < m_energySize; ++i) {
    auto E = m_eAxis[i] / SI::eV;
    auto Ip = SI::cLight / (4.0 * M_PI) * m_np[i] / Ip_units;
    auto E2Inu = pow2(E) * SI::cLight / (4.0 * M_PI) * m_nnu[i] / E2Inu_units;
    auto Qnu = m_qnu[i] / q_units;
    out << std::setprecision(5) << E << " " << Ip << " " << E2Inu << " " << Qnu << "\n";
  }

  LOGD << "Dumped spectrum to " << filename;
}

}  // namespace gnuprop
