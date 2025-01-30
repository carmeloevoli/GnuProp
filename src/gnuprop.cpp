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
       << "Omega_L: " << m_cosmology->OmegaL;
}

void GnuProp::build() {
  m_nuProductionRate = std::make_unique<gnuprop::NeutrinoProductionRate>();
  m_losses = std::make_unique<gnuprop::EnergyLosses>();
  if (m_doPhotoPion) {
    m_losses->enablePhotoPion();
  } else {
    m_losses->disablePhotoPion();
  }
  m_eAxis = simprop::utils::LogAxis(m_energyMin, m_energyMax, m_energySize);
  m_np.assign(m_energySize, 0.0);
  m_nnu.assign(m_energySize, 0.0);
  m_qnu.assign(m_energySize, 0.0);

  knownTerm.assign(m_energySize - 1, 0.0);
  diagonal.assign(m_energySize - 1, 0.0);
  upperDiagonal.assign(m_energySize - 2, 0.0);
  lowerDiagonal.assign(m_energySize - 2, 0.0);
}

// void GnuProp::evolve(double zObs) {
//   assert(zObs < m_zMax);
//   auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
//   std::reverse(zAxis.begin(), zAxis.end());
//   const auto dz = zAxis[0] - zAxis[1];

//   std::vector<double> npUp(m_energySize);

//   for (const auto& z : zAxis) {
//     // for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {  // Reverse iteration
//     // auto z = *it;

//     LOGD << std::setprecision(5) << z << "\t" << std::accumulate(m_np.begin(), m_np.end(), 0.);

//     const auto dtdz = std::abs(m_cosmology->dtdz(z));
//     const auto H = std::abs(m_cosmology->hubbleRate(z));

//     // evolve protons
//     //    #pragma omp parallel for  // Optional: Enable parallelism
//     for (size_t i = 0; i < m_energySize - 1.; ++i) {
//       const auto E = m_eAxis[i];
//       const auto Eup = m_eAxis[i + 1];
//       const auto Q = std::pow(1. + z, 3.0) * pComovingSource(E, z);

//       const auto b = E * (m_losses->beta(E, z) + H);
//       const auto bUp = Eup * (m_losses->beta(Eup, z) + H);
//       const auto dbndE = (bUp * m_np[i + 1] - b * m_np[i]) / (Eup - E);

//       auto value = m_np[i] + dz * (dtdz * Q - 3.0 / (1.0 + z) * m_np[i] + dtdz * dbndE);
//       value = std::max(value, 0.0);
//       npUp[i] = value;
//     }

//     m_np = std::move(npUp);
//   }
// }

void GnuProp::evolve(double zObs) {
  assert(zObs < m_zMax);
  auto zAxis = simprop::utils::LinAxis(zObs, m_zMax, m_zSize);
  const auto dz = zAxis[1] - zAxis[0];

  for (auto it = zAxis.rbegin(); it != zAxis.rend(); ++it) {
    auto z = *it;
    LOGD << std::setprecision(5) << z << "\t" << std::accumulate(m_np.begin(), m_np.end(), 0.);

    const auto dtdz = std::abs(m_cosmology->dtdz(z));
    const auto H = std::abs(m_cosmology->hubbleRate(z));

    // {
    //   // evolve protons
    //   std::vector<double> npUp(m_energySize, 0.);

    //   for (size_t i = 0; i < m_energySize - 1.; ++i) {
    //     const auto E = m_eAxis[i];
    //     const auto Eup = m_eAxis[i + 1];
    //     const auto Q = dtdz * std::pow(1. + z, 3.0) * Source(E, z) - 3. / (1. + z) * m_np[i];
    //     const auto b = E * (m_losses->beta(E, z) + H);
    //     const auto bUp = Eup * (m_losses->beta(Eup, z) + H);
    //     const auto dbndE = dtdz * (bUp * m_np[i + 1] - b * m_np[i]) / (Eup - E);
    //     const auto value = m_np[i] + dz * (Q + dbndE);
    //     assert(value > 0.);
    //     npUp[i] = value;
    //   }

    //   m_np = std::move(npUp);
    // }

    {
      std::vector<double> npUp(m_energySize - 1, 0.);

      for (size_t i = 0; i < m_energySize - 1; ++i) {
        const auto E = m_eAxis[i];
        const auto Eup = m_eAxis[i + 1];
        const auto dE = Eup - E;
        const auto Q = dtdz * std::pow(1. + z, 3.0) * Source(E, z) - 3. / (1. + z) * m_np[i];
        const auto b = dtdz * E * (m_losses->beta(E, z) + H);
        const auto bUp = dtdz * Eup * (m_losses->beta(Eup, z) + H);

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
    computeNuEmissivity(z);

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

double GnuProp::Source(double E, double z) const {
  const auto K = m_sourceComovingEmissivity * (m_injSlope - 2.0) / std::pow(m_energyMin, 2.);
  auto value = (E >= m_energyMin) ? std::pow(E / m_energyMin, -m_injSlope) : 0.;
  value *= (m_expCutoff > 0.0) ? std::exp(-E / m_expCutoff) : 1.;
  value *= (z <= m_zMax) ? std::pow(1.0 + z, m_evolutionIndex) : 0.;
  return K * std::max(value, 0.0);
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

void GnuProp::computeNuEmissivity(double z) {
  const auto ln_eRatio = std::log(m_eAxis[1] / m_eAxis[0]);
  for (size_t i = 0; i < m_energySize; ++i) {
    const auto Enu = m_eAxis[i];
    double value = 0.;
    for (size_t j = i; j < m_energySize; ++j) {
      value += m_np[j] * m_nuProductionRate->get(Enu, m_eAxis[j], z);
    }
    m_qnu[i] = ln_eRatio * value;
  }
}

}  // namespace gnuprop

// std::vector<double> computeNuEmissivity(double z) const;
// double nuInteractionRate(double Enu, double Ep, double z, size_t N = 5) const;

// std::vector<double> m_nnu;
// std::vector<double> m_Qnu;

// m_nnu = std::vector<double>(m_energySize, 0.0);
// m_Qnu = std::vector<double>(m_energySize, 0.0);
//
// const double units_emissivity = 1. / SI::eV / SI::m3 / SI::sec;
// out << SI::cLight / 4. / M_PI * m_nnu.at(i) / units << " ";
// out << m_Qnu.at(i) / units_emissivity << " ";

//