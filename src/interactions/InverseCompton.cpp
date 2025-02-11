#include "interactions/InverseCompton.h"

#include "simprop.h"

namespace Interactions {

double InverseCompton::sigma_com(double s) const {
  auto value = 0.;
  const auto x = s / pow2(SI::electronMassC2);
  if (x > 1.) {
    const auto b_e = (x - 1.) / (x + 1.);
    const auto be_2 = pow2(b_e);
    const auto be_3 = pow3(b_e);
    const auto c_1 = 2. * (2. + 2. * b_e - be_2 - 2. * be_3) / b_e / (1. + b_e);
    const auto c_2 = (2. - 3. * be_2 - be_3) / be_2 * std::log(x);
    value = 3. / 8. * SI::sigmaTh / x / b_e * (c_1 - c_2);
  }
  return std::max(value, 0.);
}

double InverseCompton::sigma_lab(double eElectron, double eBkg, double mu) const {
  const auto me_2 = pow2(SI::electronMassC2);
  const auto s = me_2 + 2. * eElectron * eBkg - 2. * std::sqrt(pow2(eElectron) - me_2) * eBkg * mu;

  return sigma_com(s);
}

double InverseCompton::dsigma_dE(double eElectron, double eBkg, double eGamma) const {
  const auto gamma_e = std::sqrt(1. + pow2(eElectron / SI::electronMassC2));
  const auto Gamma_e = 4. * eBkg * gamma_e / SI::electronMassC2;
  const auto E_1 = eGamma / gamma_e / SI::electronMassC2;
  const auto q = E_1 / Gamma_e / (1. - E_1);
  if (E_1 < eBkg / eElectron or E_1 > Gamma_e / (Gamma_e + 1)) return 0.;

  const auto c_1 = 2. * q * std::log(q) + (1. + 2. * q) * (1. - q);
  const auto c_2 = 0.5 * pow2(Gamma_e * q) * (1. - q) / (1. + Gamma_e * q);
  const auto value = 3. / 4. * SI::sigmaTh / eBkg / gamma_e * (c_1 + c_2);

  return value;
}

}  // namespace Interactions