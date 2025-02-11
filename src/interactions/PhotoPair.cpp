#include "interactions/PhotoPair.h"

#include <cmath>

#include "simprop.h"

namespace Interactions {

double PhotoPair::sigma_com(double s) const {
  const auto chi = s / 4. / pow2(SI::electronMassC2);
  if (chi < 1. || chi > 1e5) return 0.;

  const auto beta = sqrt(1. - 1. / chi);
  return 3. / 16. * SI::sigmaTh * (1. - pow2(beta)) *
         (2. * beta * (pow2(beta) - 2.) + (3. - pow4(beta)) * log((1. + beta) / (1. - beta)));
}

double PhotoPair::sigma_lab(double eGamma, double eBkg, double mu) const {
  const auto s = 2. * eGamma * eBkg * (1. - mu);
  return sigma_com(s);
}

// double PhotoPair::integrateXsec(double x) const {  // x = eps * E_gamma
//   const auto sMin = 4. * pow2(SI::electronMassC2);
//   const auto sMax = 4. * x;
//   auto integrand = [&](double s) { return s * sigmaInCoMFrame(s); };
//   auto value = simprop::utils::QAGIntegration<double>(integrand, sMin, sMax, 1000, 1e-4);
//   return value;
// }

double PhotoPair::dsigma_dE(double eGamma, double eBkg, double eLepton) const {
  const auto A = eBkg + eGamma;
  const auto me_c2 = SI::electronMassC2;
  const auto me2_c4 = pow2(SI::electronMassC2);
  const auto minLeptonEnergy = 0.5 * A * (1. - std::sqrt(1. - me2_c4 / eBkg / eGamma));
  const auto maxLeptonEnergy = 0.5 * A * (1. + std::sqrt(1. - me2_c4 / eBkg / eGamma));
  double value = 0.;
  if (eLepton > minLeptonEnergy and eLepton < maxLeptonEnergy) {
    auto c = 3. / 64. * SI::sigmaTh * me2_c4 / pow2(eBkg) / pow3(eGamma);
    auto a1 = 4. * A * me_c2 / ((A - eLepton) * eLepton) *
              std::log(4. * eBkg * (A - eLepton) * eLepton / me2_c4 / A);
    auto a2 = -8. * eBkg * A / me2_c4;
    auto a3 = 2. * (2. * eBkg * A / me2_c4 - 1.) * pow2(A) / ((A - eLepton) * eLepton);
    auto a4 = -(1. - me2_c4 / eBkg / A) * (pow4(A) / pow2(A - eLepton) / pow2(eLepton));
    value = c * (a1 + a2 + a3 + a4);
  }
  return value;
}

}  // namespace Interactions
