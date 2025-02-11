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

// eGamma + eBkg -> eLepton + ...

double PhotoPair::dsigma_dE(double eGamma, double eBkg, double eLepton) const {
  const auto _eGamma = eGamma / SI::electronMassC2;
  const auto _eBkg = eBkg / SI::electronMassC2;
  const auto _eLepton = eLepton / SI::electronMassC2;

  const auto A = _eBkg + _eGamma;
  const auto minLeptonEnergy = 0.5 * A * (1. - std::sqrt(1. - 1. / _eBkg / _eGamma));
  const auto maxLeptonEnergy = 0.5 * A * (1. + std::sqrt(1. - 1. / _eBkg / _eGamma));

  double value = 0.;

  if (_eLepton > minLeptonEnergy and _eLepton < maxLeptonEnergy) {
    auto c = 3. / 64. * SI::sigmaTh / pow2(_eBkg) / pow3(_eGamma);
    auto a1 =
        4. * A / ((A - _eLepton) * _eLepton) * std::log(4. * _eBkg * (A - _eLepton) * _eLepton / A);
    auto a2 = -8. * _eBkg * A;
    auto a3 = 2. * (2. * _eBkg * A - 1.) * pow2(A) / (A - _eLepton) / _eLepton;
    auto a4 = -(1. - 1. / _eBkg / A) * pow4(A) / pow2(A - _eLepton) / pow2(_eLepton);
    value = c * (a1 + a2 + a3 + a4);
  }

  return value / SI::electronMassC2;
}

}  // namespace Interactions
