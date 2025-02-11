#ifndef GNUPROP_INTERACTIONS_ICS_H
#define GNUPROP_INTERACTIONS_ICS_H

namespace Interactions {

struct InverseCompton {
  double sigma_com(double s) const;
  double sigma_lab(double eLepton, double eBkg, double mu) const;
  double dsigma_dE(double eLepton, double eBkg, double eGamma) const;
};

}  // namespace Interactions

#endif  // GNUPROP_INTERACTIONS_ICS_H