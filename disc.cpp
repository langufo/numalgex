#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "BndryProp.hpp"
#include "CylinderPDE.hpp"
#include "GSeidelSolver.hpp"
#include "JacobiSolver.hpp"
#include "Matrix.hpp"
#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"
#include "SORSolver.hpp"

/* rho / 4 epsilon = 1 */

Real
pot_axis(Real z, Real r2, Real h)
{
  z = std::abs(z);
  Real a = z - h / 2;
  Real b = z + h / 2;
  Real p = std::sqrt(r2 + a * a);
  Real q = std::sqrt(r2 + b * b);

  if (a < 0) {
    return b * q - a * p - r2 * std::log(r2 / ((p - a) * (q + b))) - 2 * z * z -
           h * h / 2;
  } else {
    return 2 * h * z * z / (p + q) + h / 2 * (p + q) - 2 * h * z -
           r2 * std::log((p + a) / (q + b));
  }
}

/*
 * Argument list:
 *  - matrix size
 *  - number of SOR iterations
 *  - number of iterations with another algorithm
 *  - the other solving algorithm (1 -> Gauss-Seidel, 2 -> Jacobi)
 *  - filenames prefix
 */
int
main(int argc, char * argv[])
{
  using std::ofstream;
  using std::scientific;
  using std::string;
  using std::vector;

  long n = std::strtol(argv[1], nullptr, 0);
  long iterFast = std::strtol(argv[2], nullptr, 0);
  long iterSlow = std::strtol(argv[3], nullptr, 0);
  long algo = std::strtol(argv[4], nullptr, 0);

  string prefix(argv[5]);
  ofstream iterFile(prefix + "iter.txt");
  ofstream axisFile(prefix + "axis.txt");
  ofstream solFile(prefix + "sol.txt");
  ofstream errFile(prefix + "err.txt");
  ofstream relFile(prefix + "rel.txt");
  scientific(iterFile);
  scientific(axisFile);
  scientific(solFile);
  scientific(errFile);
  scientific(relFile);

  vector<Real> zero(n);
  vector<Real> nonzero(n);
  vector<Real> sol(n * n);
  vector<Real> rhs(n * n);
  Matrix<Real> m(n, rhs.data());

  Real h = static_cast<Real>(1) / (n + 1);
  Real h2 = h * h;

  long a = n / 8 + 1;
  long l = n / 32 + 1;

  long a2 = a * a;

  for (long i = 0; i < n; ++i) {
    nonzero[i] =
      a2 * l * 2 * h2 / std::sqrt((i + 1) * (i + 1) + (n + 1) * (n + 1));
  }

  for (long i = 0; i < n; ++i) {
    for (long j = 0; j < n; ++j) {
      if (i + 1 > a || j + 1 > l) {
        m[i][j] = 0;
      } else {
        m[i][j] = -4;
      }
      if (i + 1 == a) {
        m[i][j] /= 2;
      }
      if (j + 1 == l) {
        m[i][j] /= 2;
      }
    }
  }

  PDE pde;
  pde.next_value = CylinderPDE::next_value;
  pde.residual = CylinderPDE::residual;
  pde.bottom = zero.data();
  pde.left = zero.data();
  pde.right = nonzero.data();
  pde.top = nonzero.data();
  pde.neum = BOTTOMBNDRY | LEFTBNDRY;
  pde.rhs = rhs.data();

  Real w = 2 / (1 + 2 * std::acos(static_cast<Real>(0)) / (n + 2));

  SORSolver solvFast(h, h, h, n, n, w);

  PDESolver * otherSolv = nullptr;

  switch (algo) {
    case 1:
      otherSolv = new GSeidelSolver(h, h, h, n, n);
      break;
    case 2:
      otherSolv = new JacobiSolver(h, h, h, n, n);
      break;
  }

  long iter[2] = { iterFast, iterSlow };
  PDESolver * solv[2] = { &solvFast, otherSolv };

  Real thrErr = 1;
  Real thrRes = 1;
  long globalIter = 1;
  m.set_first_elem(sol.data());
  for (int step = 0; step < 2; ++step) {
    for (long k = 0; k < iter[step]; ++k) {
      Real err = solv[step]->iter(sol.data(), pde);
      Real res = solv[step]->abs_res_sum(sol.data(), pde);
      iterFile << globalIter << "\t" << err << "\t" << res << "\t";

      Real max = 0;
      Real rel = 0;
      for (long j = 0; j < n; ++j) {
        Real e = std::abs(m[0][j] - pot_axis((j + 1) * h, a2 * h2, 2 * l * h));
        if (e > max) {
          max = e;
        }
        e /= std::abs(m[0][j]);
        if (e > rel) {
          rel = e;
        }
      }

      iterFile << max << "\t" << rel << "\n";

      /* print solution and errors */
      if (err < thrErr || res < thrRes) {
        if (err < thrErr) {
          thrErr /= 10;
        }
        if (res < thrRes) {
          thrRes /= 10;
        }

        solFile << "#\t" << globalIter << "\t" << err << "\t" << res << "\n";
        axisFile << "#\t" << globalIter << "\t" << err << "\t" << res << "\n";
        errFile << "#\t" << globalIter << "\t" << err << "\t" << res << "\n";
        relFile << "#\t" << globalIter << "\t" << err << "\t" << res << "\n";
        for (long j = 0; j < n; ++j) {
          Real y = (j + 1) * h;

          axisFile << y << "\t" << m[0][j] << "\n";

          for (long i = 0; i < n; ++i) {
            solFile << (i + 1) * h << "\t" << y << "\t" << m[i][j] << "\n";
          }
          solFile << "\n";

          Real e = std::abs(m[0][j] - pot_axis(y, a2 * h2, 2 * l * h));
          errFile << y << "\t" << e << "\n";
          e /= std::abs(m[0][j]);
          relFile << y << "\t" << e << "\n";
        }
        solFile << "\n";
        axisFile << "\n\n";
        errFile << "\n\n";
        relFile << "\n\n";
      }

      ++globalIter;
    }
    iterFile << "\n";
  }

  delete otherSolv;
}
