#ifndef GNUPROP_INTERACTIONS_BREITWHEELER_H
#define GNUPROP_INTERACTIONS_BREITWHEELER_H

namespace PhotonPairProduction {

struct BreitWheeler {
  double sigmaInCoMFrame(double s) const;
  double sigma(double eGamma, double eBkg, double mu) const;

  double integrateXsec(double x) const;

  double dsigmadE(double eGamma, double eLepton, double eBkg) const;
};

}  // namespace PhotonPairProduction

#endif  // GNUPROP_INTERACTIONS_BREITWHEELER_H