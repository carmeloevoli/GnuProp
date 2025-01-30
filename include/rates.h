#ifndef GNUPROP_RATES_H
#define GNUPROP_RATES_H

#include <string>
#include <vector>

namespace gnuprop {

class NeutrinoProductionRate {
 public:
  NeutrinoProductionRate();

  // Compute beta
  double get(double E_nu, double E_p, double z = 0) const;

 private:
  const std::string m_filename = "data/gnuprop_nu_production_test.bin";

  const double m_zMin = 0;
  const double m_zMax = 10;
  const size_t m_zSize = 51;

  const double m_lgXMin = -4;
  const double m_lgXMax = 0;
  const size_t m_xSize = 200;

  const double m_lgEnergyMin = 19;
  const double m_lgEnergyMax = 22;
  const size_t m_energySize = 400;

  std::vector<double> m_z;
  std::vector<double> m_lgx;
  std::vector<double> m_lgE;

  std::vector<std::vector<double>> m_rate;
};

class GammaAbsorptionRate {
 public:
  GammaAbsorptionRate();

  // Compute beta
  double get(double E_gamma, double z = 0) const;

 private:
  const std::string m_filename = "data/gnuprop_gamma_absorption_rate.bin";

  const double m_zMin = 0;
  const double m_zMax = 10;
  const size_t m_zSize = 101;

  const double m_lgEnergyMin = 13;
  const double m_lgEnergyMax = 22;
  const size_t m_energySize = 3000;

  std::vector<double> m_z;
  std::vector<double> m_lgE;

  std::vector<std::vector<double>> m_rate;
};

}  // namespace gnuprop

#endif  // GNUPROP_RATES_H