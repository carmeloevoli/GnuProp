#ifndef GNUPROP_LOSSES_H
#define GNUPROP_LOSSES_H

#include <string>
#include <vector>

namespace gnuprop {

class EnergyLosses {
 public:
  EnergyLosses();

  // Compute beta
  double beta(double E, double z = 0) const;

  // PhotoPion control
  void enablePhotoPion() { doPhotoPion = true; }
  void disablePhotoPion() { doPhotoPion = false; }
  bool isPhotoPionEnabled() const { return doPhotoPion; }

 private:
  bool doPhotoPion = false;
  const size_t m_energySize = 10000;
  const double m_lgEnergyMin = 16;
  const double m_lgEnergyMax = 25;
  const std::string m_filename_pair = "data/gnuprop_proton_pair_losses_16_25_1e4.bin";
  const std::string m_filename_photopi = "data/gnuprop_proton_photopion_losses_16_25_1e4.bin";

  std::vector<double> m_lgE;
  std::vector<double> m_beta_pair;
  std::vector<double> m_beta_photopion;
};

}  // namespace gnuprop

#endif  // GNUPROP_LOSSES_H
