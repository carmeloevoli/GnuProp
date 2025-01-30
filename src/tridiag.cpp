#include "tridiag.h"

void solveTridiagonalSystem(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c,
                            std::vector<double>& r, std::vector<double>& u, int n) {
  int j = 0;
  double bet = 0.0;
  std::vector<double> gam(n, 0.);  // double gam[n];
  //	One vector of workspace, gam, is needed.
  if (b[0] == 0.0) std::cerr << "Error 1 in tridag: the first diagonal term is 0!! " << std::endl;
  // If this happens, then you should rewrite your equations as a set of order N-1, with u1
  // trivially eliminated.
  bet = b[0];
  u[0] = r[0] / bet;
  for (j = 1; j < n; j++) {  // Decomposition and forward substitution.
    // double* gm = gam+j;
    //(*gm) = c[j-1]/bet;
    gam[j] = c[j - 1] / bet;
    // bet = b[j] - a[j]*(*gm);
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0) {
      std::cout << "j = 0 " << " --> diagonal term b[0] = " << b[0]
                << " off diagonal term a[0] = " << a[0] << " c[0] = " << c[0] << " u[0] = " << u[0]
                << " bet = b[0] " << std::endl;
      std::cout << "j = " << j << " --> diagonal term b[j] = " << b[j]
                << " off diagonal term a[j] = " << a[j] << " gam[j] = " << gam[j]
                << " bet = b[j] - a[j]*c[j-1]/bet " << bet << std::endl;
      std::cerr << "Error 2 in tridag: bet = 0!" << std::endl;
    }
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }
  for (j = (n - 2); j >= 0; j--) u[j] -= gam[j + 1] * u[j + 1];  // Backsubstitution.
  return;
}

/* plain gauss elimination, only not bothering with the zeroes
 *
 *       diag[0]  abovediag[0]             0   .....
 *  belowdiag[0]       diag[1]  abovediag[1]   .....
 *             0  belowdiag[1]       diag[2]
 *             0             0  belowdiag[2]   .....
 */
int solve_tridiag_nonsym(const std::vector<double>& diag, const std::vector<double>& abovediag,
                         const std::vector<double>& belowdiag, const std::vector<double>& rhs,
                         std::vector<double>& x, size_t N) {
  int status = GSL_SUCCESS;
  std::vector<double> alpha(N);
  std::vector<double> z(N);

  size_t i, j;

  /* Bidiagonalization (eliminating belowdiag)
   & rhs update
   diag' = alpha
   rhs' = z
   */
  alpha[0] = diag.at(0);
  z[0] = rhs.at(0);

  if (alpha[0] == 0) {
    status = GSL_EZERODIV;
  }

  for (i = 1; i < N; i++) {
    const double t = belowdiag.at(i - 1) / alpha[i - 1];
    alpha[i] = diag.at(i) - t * abovediag.at(i - 1);
    z[i] = rhs.at(i) - t * z[i - 1];
    if (alpha[i] == 0) {
      status = GSL_EZERODIV;
    }
  }

  /* backsubstitution */
  x.at(N - 1) = z[N - 1] / alpha[N - 1];
  if (N >= 2) {
    for (i = N - 2, j = 0; j <= N - 2; j++, i--) {
      x.at(i) = (z[i] - abovediag.at(i) * x.at(i + 1)) / alpha[i];
    }
  }

  if (status == GSL_EZERODIV) {
    std::cout << "Error : matrix must be positive definite!" << "\n";
  }

  // delete alpha;
  // delete z;

  return status;
}

int gsl_linalg_solve_tridiag(const std::vector<double>& diag, const std::vector<double>& abovediag,
                             const std::vector<double>& belowdiag, const std::vector<double>& rhs,
                             std::vector<double>& solution) {
  if (diag.size() != rhs.size()) {
    std::cout << "size of diag must match rhs" << std::endl;
    exit(GSL_EBADLEN);
  } else if (abovediag.size() != rhs.size() - 1) {
    std::cout << "size of abovediag must match rhs-1" << std::endl;
    exit(GSL_EBADLEN);
  } else if (belowdiag.size() != rhs.size() - 1) {
    std::cout << "size of belowdiag must match rhs-1" << std::endl;
    exit(GSL_EBADLEN);
  } else if (solution.size() != rhs.size()) {
    std::cout << "size of solution must match rhs" << std::endl;
    exit(GSL_EBADLEN);
  } else {
    return solve_tridiag_nonsym(diag, abovediag, belowdiag, rhs, solution, diag.size());
  }
  return 0;
}
