#include "interactions/BreitWheeler.h"

#include <cmath>

#include "simprop.h"

namespace BreitWheeler {

double OpticalDepth::sigmaInCoMFrame(const double &s) const {
  const auto chi = s / 4. / pow2(SI::electronMassC2);

  if (chi < 1. || chi > 1e5) return 0.;

  const auto beta = sqrt(1. - 1. / chi);
  return 3. / 16. * SI::sigmaTh * (1. - pow2(beta)) *
         (2. * beta * (pow2(beta) - 2.) + (3. - pow4(beta)) * log((1. + beta) / (1. - beta)));
}

double OpticalDepth::sigma(const double &eGamma, const double &eBkg, const double &mu) const {
  const auto s = 2. * eGamma * eBkg * (1. - mu);
  return sigmaInCoMFrame(s);
}

double OpticalDepth::integrateXsec(double x) const {  // x = eps * E_gamma
  const auto sMin = 4. * pow2(SI::electronMassC2);
  const auto sMax = 4. * x;
  auto integrand = [&](double s) { return s * sigmaInCoMFrame(s); };
  auto value = simprop::utils::QAGIntegration<double>(integrand, sMin, sMax, 1000, 1e-4);
  return value;
}

}  // namespace BreitWheeler
