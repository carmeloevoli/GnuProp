#ifndef GNUPROP_INTERACTIONS_PHOTOPAIR_H
#define GNUPROP_INTERACTIONS_PHOTOPAIR_H

namespace Interactions {

struct PhotoPair {
  double sigma_com(double s) const;
  double sigma_lab(double eGamma, double eBkg, double mu) const;
  double dsigma_dE(double eGamma, double eBkg, double eLepton) const;
};

}  // namespace Interactions

#endif  // GNUPROP_INTERACTIONS_PHOTOPAIR_H