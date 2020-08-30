#include <iostream>
#include <math.h>

#include "GSeidelCyl.hpp"
#include "JacobiCyl.hpp"
#include "Matrix.hpp"
#include "PDESolver.hpp"
#include "SORCyl.hpp"

const int n = 512;
int hHalfN = 16;
int rIntN = 32;
int rExtN = 64;
double h = 1.0 / n;
double height = hHalfN * 2 * h;
double intRadius = rIntN * h;
double extRadius = rExtN * h;
double density = 1;

double
ax_pot(double z)
{
  z = fabs(z);
  double rInt2 = intRadius * intRadius;
  double rExt2 = extRadius * extRadius;
  double a = (z - 0.5 * height);
  double b = (z + 0.5 * height);
  double p = sqrt(rInt2 + a * a);
  double q = sqrt(rInt2 + b * b);
  double s = sqrt(rExt2 + a * a);
  double t = sqrt(rExt2 + b * b);
  return -0.25 * density *
         (rExt2 * log((a + s) / (b + t)) - rInt2 * log((a + p) / (b + q)) +
          (rExt2 - rInt2) * (a / (s + p) - b / (t + q)));
}

int
main()
{
  double tol = 6e-5;

  Matrix<double> m(n);

  double rhs[n * n];
  m.set_first_elem(rhs);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i + 1 < rIntN || i + 1 > rExtN || j + 1 > hHalfN)
        m[i][j] = 0;
      else
        m[i][j] = -density;

      if (i + 1 == rIntN)
        m[i][j] *= 0.5;
      if (i + 1 == rExtN)
        m[i][j] *= 0.5;
      if (j + 1 == hHalfN)
        m[i][j] *= 0.5;
    }
  }

  double b[n], t[n], l[n], r[n];
  for (int i = 0; i < n; ++i) {
    b[i] = 0;
    l[i] = 0;

    long long r2 = (i + 1) * (i + 1) + (n + 1) * (n + 1);
    double a =
      0.25 * density * height * (rExtN * rExtN - rIntN * rIntN) * h / sqrt(r2);
    t[i] = a;
    r[i] = a;
  }

  PDESolver::BndryLayout neum = PDESolver::BOTTOMBNDRY | PDESolver::LEFTBNDRY;
  SORCyl sorSolv(rhs, neum, h, h, h, n, n, b, t, l, r,
                 2 * n / (n + 2 * acos(0)));
  sorSolv.rev_solv_direc(true, true);

  double s[n * n] = {};
  m.set_first_elem(s);
  int iter = 0;
  double eps;
  do {
    eps = sorSolv.iter(s, true);
    ++iter;
    if (iter % 128 == 0) {
      std::cerr << eps << "\n";
    }
  } while (eps > tol);

  using std::cout;
  cout << "# eps = " << eps << "\n";
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cout << (i + 1) * h << " " << (j + 1) * h << " " << m[i][j];

      if (i == 0) {
        cout << " " << m[i][j] << " " << ax_pot((j + 1) * h) << " "
             << m[i][j] - ax_pot((j + 1) * h);
      }

      cout << "\n";
    }
    cout << "\n";
  }
}
