#include <memory>

#include "gnuprop.h"
#include "simprop.h"

class TestProtons {
 public:
  static std::unique_ptr<gnuprop::GnuProp> testEvolutionIndex(double m) {
    auto g = std::make_unique<gnuprop::GnuProp>(std::make_unique<simprop::cosmo::Cosmology>());
    g->setEvolutionIndex(m);
    g->build();
    g->evolve(0.);
    return g;
  }

  static std::unique_ptr<gnuprop::GnuProp> testInjectionSlope(double gamma) {
    auto g = std::make_unique<gnuprop::GnuProp>(std::make_unique<simprop::cosmo::Cosmology>());
    g->setInjectionSlope(gamma);
    g->build();
    g->evolve(0.);
    return g;
  }

  static std::unique_ptr<gnuprop::GnuProp> testzmax(double zMax) {
    auto g = std::make_unique<gnuprop::GnuProp>(std::make_unique<simprop::cosmo::Cosmology>());
    g->setRedshiftMax(zMax);
    g->build();
    g->evolve(0.);
    return g;
  }
};

int main() {
  try {
    simprop::utils::startup_information();
    // auto base = TestProtons::testEvolutionIndex(0.);
    // base->dump("gnuprop_protons_test_base.txt");
    // auto m3 = TestProtons::testEvolutionIndex(3.);
    // m3->dump("gnuprop_protons_test_m3.txt");
    // auto m5 = TestProtons::testEvolutionIndex(5.);
    // m5->dump("gnuprop_protons_test_m5.txt");
    // auto m3neg = TestProtons::testEvolutionIndex(-3.);
    // m3neg->dump("gnuprop_protons_test_m-3.txt");
    // auto m5neg = TestProtons::testEvolutionIndex(-5.);
    // m5neg->dump("gnuprop_protons_test_m-5.txt");

    // auto g21 = TestProtons::testInjectionSlope(2.1);
    // g21->dump("gnuprop_protons_test_gamma21.txt");
    // auto g23 = TestProtons::testInjectionSlope(2.3);
    // g23->dump("gnuprop_protons_test_gamma23.txt");
    // auto g27 = TestProtons::testInjectionSlope(2.7);
    // g27->dump("gnuprop_protons_test_gamma27.txt");
    // auto g30 = TestProtons::testInjectionSlope(3.0);
    // g30->dump("gnuprop_protons_test_gamma30.txt");

    auto g1 = TestProtons::testzmax(1.);
    g1->dump("gnuprop_protons_test_z1.txt");
    auto g2 = TestProtons::testzmax(2.);
    g2->dump("gnuprop_protons_test_z2.txt");
    auto g3 = TestProtons::testzmax(3.);
    g3->dump("gnuprop_protons_test_z3.txt");
    auto g4 = TestProtons::testzmax(4.);
    g4->dump("gnuprop_protons_test_z4.txt");
    auto g5 = TestProtons::testzmax(5.);
    g5->dump("gnuprop_protons_test_z5.txt");
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what();
  }
}