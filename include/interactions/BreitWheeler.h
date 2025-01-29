#ifndef GNUPROP_INTERACTIONS_BREITWHEELER_H
#define GNUPROP_INTERACTIONS_BREITWHEELER_H

namespace BreitWheeler {

struct OpticalDepth {
  double sigmaInCoMFrame(const double &s) const;
  double sigma(const double &eGamma, const double &eBkg, const double &mu) const;

  double integrateXsec(double x) const;
};

}  // namespace BreitWheeler

#endif  // GNUPROP_INTERACTIONS_BREITWHEELER_H