#include "interactions/BreitWheeler.h"
#include "simprop.h"

class OpticalDepth {
 public:
  OpticalDepth(std::unique_ptr<simprop::photonfields::PhotonField> phField)
      : m_phField(std::move(phField)) {}
  virtual ~OpticalDepth() = default;

  double get(double E_0, double z) const {
    auto integrateOverField = [&](double E_0, double z) {
      const auto epsMin = m_phField->getMinPhotonEnergy();
      const auto epsMax = m_phField->getMaxPhotonEnergy();

      auto integrand = [&](double lnEpsilon) {
        const auto epsilon = std::exp(lnEpsilon);
        const auto x = E_0 * epsilon * (1. + z);
        if (x > pow2(SI::electronMassC2))
          return m_phField->density(epsilon, z) / epsilon * m_bw.intSSigma(x);
        else
          return 0.;
      };

      const auto a = std::log(epsMin);
      const auto b = std::log(epsMax);
      const size_t N = 10000;
      auto value = simprop::utils::QAGIntegration<double>(integrand, a, b, N, 1e-3);
      return value;
    };

    auto integrand = [&](double z) {
      return 1. / pow2(1. + z) * m_cosmology.dtdz(z) * integrateOverField(E_0, z);
    };
    auto value = simprop::utils::QAGIntegration<double>(integrand, 0., z, 10000, 1e-3);
    return SI::cLight / pow2(E_0) / 8. * value;
  }

 private:
  std::unique_ptr<simprop::photonfields::PhotonField> m_phField;
  PhotonPairProduction::BreitWheeler m_bw;
  simprop::cosmo::Planck2018 m_cosmology;
};

void dumpOpticalDepth(std::unique_ptr<simprop::photonfields::PhotonField> phField,
                      std::string filename) {
  const auto redshifts = std::vector<double>({0.1, 0.5, 1.0, 3.0, 5.0, 10.});
  const auto eGamma = simprop::utils::LogAxis<double>(1e-2 * SI::TeV, 1e6 * SI::TeV, 100);
  OpticalDepth opticalDepth(std::move(phField));

  std::ofstream out(filename);
  out << "# energy [TeV] - tau\n";
  out << std::scientific;
  for (const auto& E : eGamma) {
    out << E / SI::TeV << "\t";
    for (const auto& z : redshifts) out << opticalDepth.get(E, z) << "\t";
    out << "\n";
  }
}

void dumpPhotonField(std::unique_ptr<simprop::photonfields::PhotonField> phField,
                     std::string filename) {
  const auto redshifts = std::vector<double>({0., 0.1, 0.5, 1.0, 3.0, 5.0, 10.});
  const auto epsAxis = simprop::utils::LogAxis<double>(1e-5 * SI::eV, 1e5 * SI::eV, 1000);
  const auto units = 1. / SI::cm3;

  std::ofstream out(filename);
  out << "# energy [eV] - tau\n";
  out << std::scientific;
  for (const auto& eps : epsAxis) {
    out << eps / SI::eV << "\t";
    for (const auto& z : redshifts)
      out << eps * phField->density(eps, z) / pow3(1. + z) / units << "\t";
    out << "\n";
  }
}

int main() {
  try {
    simprop::utils::startup_information();

    simprop::utils::Timer timer("main timer");

    dumpPhotonField(std::make_unique<simprop::photonfields::CMB>(),
                    "output/gnuprop_phfield_cmb.txt");
    dumpPhotonField(std::make_unique<simprop::photonfields::Saldana2021PhotonField>(),
                    "output/gnuprop_phfield_ebl.txt");
    dumpOpticalDepth(std::make_unique<simprop::photonfields::CMB>(),
                     "output/gnuprop_optical_depth_cmb.txt");
    dumpOpticalDepth(std::make_unique<simprop::photonfields::Saldana2021PhotonField>(),
                     "output/gnuprop_optical_depth_ebl.txt");

  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}