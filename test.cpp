#include <iostream>
#include <math.h>

#include "GSeidelCyl.hpp"
#include "JacobiCyl.hpp"
#include "Matrix.hpp"
#include "PDESolver.hpp"
#include "SORCyl.hpp"

#define SORONLY
#define BNDRYCORR

const double l = 1;
const int n = 128;

const double h = l / (n + 1);

const double rInt = h + n / 4 * h;
const double rExt = h + n / 2 * h;

const double tol = 1e-7;

void
zero(double *s);

void
test(PDESolver &solv, double *s, const double *t);

double
pot_hollow_cyl(double r)
{
  if (r <= rInt) {
    return 0;
  } else if (r <= rExt) {
    return (r * r - rInt * rInt) / 2 + rInt * rInt * log(rInt / r);
  } else {
    return (rExt * rExt - rInt * rInt) * log(r / rExt) + pot_hollow_cyl(rExt);
  }
}

double
rhs(double r, double z)
{
  if (r < rInt) {
    return 0;
#ifdef BNDRYCORR
  } else if (r == rInt) {
    return 1;
#endif
  } else if (r < rExt) {
    return 2;
#ifdef BNDRYCORR
  } else if (r == rExt) {
    return 1;
#endif
  } else {
    return 0;
  }
}

int
main()
{
  const double pi = 2 * acos(0);

  std::cout << std::scientific;

  double w = 2 / (1 + pi / (n + 2));
  double b[n], t[n], l[n], r[n];

  double truth[n * n];

  std::cout << "# Actual solution\n";
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double a = pot_hollow_cyl((i + 1) * h);
      truth[i * n + j] = a;
      std::cout << (i + 1) * h << " " << (j + 1) * h << " " << a << "\n";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  const char *intro[2] = { "Dirichlet only", "Dirichlet and Neumann" };

  double rhsm[n * n] = {};
  Matrix<double> m(n, rhsm);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m[i][j] = rhs((i + 1) * h, (j + 1) * h);
    }
  }

  PDESolver::BndryLayout neum = 0;
  for (int i = 0; i < n; ++i) {
    b[i] = pot_hollow_cyl(h * (i + 1));
    t[i] = pot_hollow_cyl(h * (i + 1));
    l[i] = pot_hollow_cyl(0);
    r[i] = pot_hollow_cyl(h * (n + 1));
  }

  double sol[n * n];

  for (const char *intro : intro) {
    std::cout << "# " << intro << "\n\n";

    JacobiCyl jSolv(rhsm, neum, h, h, h, n, n, b, t, l, r);
    GSeidelCyl gsSolv(rhsm, neum, h, h, h, n, n, b, t, l, r);
    SORCyl sorSolv(rhsm, neum, h, h, h, n, n, b, t, l, r, w);

    gsSolv.rev_solv_direc(true, false);

#ifndef SORONLY
    std::cout << "# Jacobi\n";
    zero(sol);
    //     test(jSolv, sol, truth);

    std::cout << "# Gauss-Seidel\n";
    zero(sol);
    test(gsSolv, sol, truth);
#endif

    std::cout << "# SOR\n";
    zero(sol);
    test(sorSolv, sol, truth);

    neum = PDESolver::BOTTOMBNDRY | PDESolver::TOPBNDRY | PDESolver::LEFTBNDRY;
    for (int i = 0; i < n; ++i) {
      b[i] = 0;
      t[i] = 0;
      l[i] = 0;
    }
  }
}

void
zero(double *s)
{
  for (int i = 0; i < n * n; ++i) {
    s[i] = 0;
  }
}

void
test(PDESolver &solv, double *s, const double *t)
{
  double h = l / (n + 1);

  double eps;
  int nIter = 0;
  do {
    eps = solv.iter(s, true);
    ++nIter;
    if (nIter % 128 == 0) {
      std::cerr << "eps = " << eps << "\n";
    }
  } while (eps > tol);

  std::cout << "# residual " << eps << "\n";
  std::cout << "# " << nIter << " iterations\n";
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      int k = i * n + j;
      std::cout << (i + 1) * h << " " << (j + 1) * h << " " << s[k] << " "
                << s[k] - t[k] << " " << (s[k] / t[k] - 1) << "\n";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}
