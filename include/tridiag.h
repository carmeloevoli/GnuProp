#ifndef TRIDIAG_H_
#define TRIDIAG_H_

#include <iostream>
#include <vector>

#define GSL_SUCCESS 0
#define GSL_ENOMEM 1
#define GSL_EZERODIV 2
#define GSL_EBADLEN 3

int gsl_linalg_solve_tridiag(const std::vector<double>& diag, const std::vector<double>& abovediag,
                             const std::vector<double>& belowdiag, const std::vector<double>& rhs,
                             std::vector<double>& solution);

int solve_tridiag_nonsym(const std::vector<double>& diag, const std::vector<double>& abovediag,
                         const std::vector<double>& belowdiag, const std::vector<double>& rhs,
                         std::vector<double>& x, size_t N);

void solveTridiagonalSystem(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c,
                            std::vector<double>& r, std::vector<double>& u, int n);

#endif /* TRIDIAG_H_ */
