#ifndef GNUPROP_INTERACTIONS_ICS_H
#define GNUPROP_INTERACTIONS_ICS_H

namespace Interactions {

struct InverseCompton {
  double sigmaInCoMFrame(double s) const;
  double sigma(double eGamma, double eBkg, double mu);
  double dsigmadE(double eGamma, double eLepton, double eBkg);
};

}  // namespace Interactions

#endif  // GNUPROP_INTERACTIONS_ICS_H