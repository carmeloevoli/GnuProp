#ifndef GNUPROP_INTERACTIONS_GAMMAPAIRPRODUCTION_H
#define GNUPROP_INTERACTIONS_GAMMAPAIRPRODUCTION_H

namespace Interactions {

struct GammaPairProduction {
  double sigmaInCoMFrame(double s) const;
  double sigma(double eGamma, double eBkg, double mu) const;
  double integrateXsec(double x) const;
  double dsigmadE(double eGamma, double eLepton, double eBkg) const;
};

}  // namespace Interactions

#endif  // GNUPROP_INTERACTIONS_GAMMAPAIRPRODUCTION_H