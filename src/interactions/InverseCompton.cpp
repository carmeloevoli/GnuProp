#include "interactions/InverseCompton.h"

#include "simprop.h"

namespace InverseCompton {

double sigma_of_s(double s) {
  const auto me_2 = pow2(SI::electronMassC2);
  if (s < me_2) return 0.;
  const auto b_e = (s - me_2) / (s + me_2);
  const auto be_2 = pow2(b_e);
  const auto be_3 = pow3(b_e);
  const auto c_0 = (1 - b_e) / (1 + b_e) / b_e;
  const auto c_1 = 2. * (2. + 2. * b_e - be_2 - 2. * be_3) / b_e / (1. + b_e);
  const auto c_2 = (2. - 3. * be_2 - be_3) / be_2 * std::log((1 + b_e) / (1. - b_e));
  const auto value = 3. / 8. * SI::sigmaTh * c_0 * (c_1 - c_2);
  return value;
}

double sigma(double eElectron, double eBkg, double mu) {
  const auto me_2 = pow2(SI::electronMassC2);
  const auto s = me_2 + 2. * eElectron * eBkg - 2. * std::sqrt(pow2(eElectron) - me_2) * eBkg * mu;
  return sigma_of_s(s);
}

double dsigmadE(double eElectron, double eBkg, double eGammaPrime) {
  const auto gamma_e = std::sqrt(1. + pow2(eElectron) / pow2(SI::electronMassC2));
  const auto Gamma_e = 4. * eBkg * gamma_e / SI::electronMassC2;
  const auto E_1 = eGammaPrime / gamma_e / SI::electronMassC2;
  const auto q = E_1 / Gamma_e / (1. - E_1);
  if (E_1 < eBkg / eElectron or E_1 > Gamma_e / (Gamma_e + 1)) return 0.;
  const auto c_1 = 2. * q * std::log(q) + (1. + 2. * q) * (1. - q);
  const auto c_2 = 0.5 * (Gamma_e * pow2(q)) * (1. - q) / (1. + Gamma_e * q);
  const auto value = 3. / 4. * SI::sigmaTh / eBkg / gamma_e * (c_1 + c_2);
  return value;
}

}  // namespace InverseCompton