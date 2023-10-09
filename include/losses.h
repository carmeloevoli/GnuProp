#ifndef BENIAMINO_LOSSES_H
#define BENIAMINO_LOSSES_H

#include <string>
#include <vector>

namespace beniamino {

class LossesTable {
 public:
  LossesTable() {}
  bool loadTable(const std::string &filename);
  double beta(double E) const;
  double dbdE(double E) const;

 private:
  std::vector<double> logE_;
  std::vector<double> logbeta_;
  std::vector<double> logdbdE_;
};

}  // namespace beniamino

#endif