#ifndef GNUPROP_INTERACTIONS_ICS_H
#define GNUPROP_INTERACTIONS_ICS_H

namespace InverseCompton {

double sigmaInCoMFrame(double s);
double sigma(double eGamma, double eBkg, double mu);

double dsigmadE(double eGamma, double eLepton, double eBkg);

}  // namespace InverseCompton

#endif  // GNUPROP_INTERACTIONS_ICS_H