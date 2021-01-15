#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

Real
pot_axis(Real z, Real a)
{
  if (z == 0) {
    return 0;
  }

  bool neg = z < 0;
  Real z2 = z * z;
  Real a2 = a * a;
  Real r2 = a2 + z2;
  Real r = std::sqrt(r2);
  z = std::abs(z);
  Real v;
  if (z > a) {
    v = (r2 * r / z - z2) / 3 - a2 / 2;
  } else {
    v = (r2 * r - a2 * a) / (3 * z) - z2 / 2;
  }
  return neg ? -v : v;
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
  using std::flush;
  using std::ofstream;
  using std::scientific;
  using std::string;
  using std::vector;

  long n = std::strtol(argv[1], nullptr, 0);
  long iterFast = std::strtol(argv[2], nullptr, 0);
  long iterSlow = std::strtol(argv[3], nullptr, 0);
  long algo = std::strtol(argv[4], nullptr, 0);

  if (n % 2 == 0) {
    std::cerr << "Dimension must be odd.\n";
    return 1;
  }

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
  vector<Real> r(n);
  vector<Real> t(n);
  vector<Real> sol(n * n);
  vector<Real> rhs(n * n);
  vector<Real> an(n);

  Real h = static_cast<Real>(1) / (n + 1);

  long a = n / 8 + 1;

  for (int i = 0; i < n; ++i) {
    an[i] = pot_axis((i + 1) * h, a * h);
    long d2 = (i + 1) * (i + 1) + (n + 1) * (n + 1);
    Real coeff = a * a * a * a * h * h / 8;
    r[i] = coeff * (i + 1) / (d2 * std::sqrt(d2));
    t[i] = coeff * (n + 1) / (d2 * std::sqrt(d2));
  }

  Matrix<Real> m(n, rhs.data());
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      long r = i + 1;
      long z = j + 1;
      if (r * r + z * z <= a * a) {
        m[i][j] = -1;
      } else {
        m[i][j] = 0;
      }
    }
  }

  PDE pde;
  pde.next_value = CylinderPDE::next_value;
  pde.residual = CylinderPDE::residual;
  pde.bottom = zero.data();
  pde.left = zero.data();
  pde.right = r.data();
  pde.top = t.data();
  pde.neum = LEFTBNDRY;
  pde.rhs = rhs.data();

  Real w = 2 / (1 + 2 * std::sqrt(4.45) / (n + 2));
  SORSolver solvFast(h, h, h, n, n, w);

  PDESolver * solvSlow = nullptr;
  switch (algo) {
    case 1:
      solvSlow = new GSeidelSolver(h, h, h, n, n);
      break;
    case 2:
      solvSlow = new JacobiSolver(h, h, h, n, n);
      break;
  }

  long iter[2] = { iterFast, iterSlow };
  PDESolver * solv[2] = { &solvFast, solvSlow };

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
        Real e = std::abs(m[0][j] - an[j]);
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
      if (err != 0 && (err < thrErr || res < thrRes)) {
        if (err < thrErr) {
          thrErr /= 10;
        }
        if (res < thrRes) {
          thrRes /= 10;
        }

        solFile << "# " << globalIter << " " << err << " " << res << "\n";
        axisFile << "# " << globalIter << " " << err << " " << res << "\n";
        errFile << "# " << globalIter << " " << err << " " << res << "\n";
        relFile << "# " << globalIter << " " << err << " " << res << "\n";

        for (long j = 0; j < n; ++j) {
          Real y = (j + 1) * h;

          axisFile << y << "\t" << m[0][j] << "\n";

          solFile << "\n";
          for (long i = 0; i < n; ++i) {
            solFile << (i + 1) * h << "\t" << y << "\t" << m[i][j] << "\n";
          }

          Real e = std::abs(m[0][j] - an[j]);
          errFile << y << "\t" << e << "\n";
          e /= std::abs(m[0][j]);
          relFile << y << "\t" << e << "\n";
        }

        solFile << "\n\n";
        axisFile << "\n\n";
        errFile << "\n\n";
        relFile << "\n\n";

        flush(solFile);
        flush(axisFile);
        flush(errFile);
        flush(relFile);
      }

      ++globalIter;
    }
    iterFile << "\n";
  }

  delete solvSlow;
}
