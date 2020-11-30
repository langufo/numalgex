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

Real
pot_sphere(Real r2, Real a2)
{
  if (r2 <= a2) {
    return (-r2 + a2 * 3) / 2;
  } else {
    return a2 * std::sqrt(a2 / r2);
  }
}

/*
 * Argument list:
 *  - matrix size
 *  - number of iterations
 *  - solving algorithm: 0 -> SOR, 1 -> Gauss-Seidel, 2 -> Jacobi
 *  - output files name prefix
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
  long iter = std::strtol(argv[2], nullptr, 0);
  long algo = std::strtol(argv[3], nullptr, 0);

  string prefix(argv[4]);
  ofstream iterFile(prefix + "iter.txt");
  ofstream axisFile(prefix + "axis.txt");
  ofstream solFile(prefix + "sol.txt");
  ofstream errFile(prefix + "err.txt");
  ofstream relFile(prefix + "rel.txt");
  scientific(iterFile);
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
  long a = n / 2 + 1;
  long a2 = a * a;

  for (int i = 0; i < n; ++i) {
    nonzero[i] =
      pot_sphere(((i + 1) * (i + 1) + (n + 1) * (n + 1)) * h2, a2 * h2);
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if ((i + 1) * (i + 1) + (j + 1) * (j + 1) <= a * a) {
        m[i][j] = -3;
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
  pde.right = nonzero.data();
  pde.top = nonzero.data();
  pde.neum = BOTTOMBNDRY | LEFTBNDRY;
  pde.rhs = rhs.data();

  PDESolver * solver = nullptr;

  Real w = 2 / (1 + 2 * std::acos(static_cast<Real>(0)) / (n + 2));

  switch (algo) {
    case 0:
      solver = new SORSolver(h, h, h, n, n, w);
      break;
    case 1:
      solver = new GSeidelSolver(h, h, h, n, n);
      break;
    case 2:
      solver = new JacobiSolver(h, h, h, n, n);
      break;
  }

  Real thrErr = 1;
  Real thrRes = 1;
  m.set_first_elem(sol.data());
  for (int k = 0; k < iter; ++k) {
    Real err = solver->iter(sol.data(), pde);
    Real res = solver->abs_res_sum(sol.data(), pde);
    iterFile << k + 1 << "\t" << err << "\t" << res << "\t";

    Real max = 0;
    Real rel = 0;
    for (long i = 0; i < n; ++i) {
      for (long j = 0; j < n; ++j) {
        Real r2 = ((i + 1) * (i + 1) + (j + 1) * (j + 1)) * h2;
        Real e = std::abs(m[i][j] - pot_sphere(r2, a2 * h2));
        if (e > max) {
          max = e;
        }
        e /= std::abs(m[i][j]);
        if (e > rel) {
          rel = e;
        }
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

      solFile << "# " << k + 1 << " " << err << " " << res << "\n";
      axisFile << "# " << k + 1 << " " << err << " " << res << "\n";
      errFile << "# " << k + 1 << " " << err << " " << res << "\n";
      relFile << "# " << k + 1 << " " << err << " " << res << "\n";

      for (long j = 0; j < n; ++j) {
        solFile << "\n";
        errFile << "\n";
        relFile << "\n";

        Real y = (j + 1) * h;
        axisFile << y << "\t"
                 << std::abs(m[0][j] - pot_sphere(h * h + y * y, a2 * h2))
                 << "\n";
        for (long i = 0; i < n; ++i) {
          Real x = (i + 1) * h;

          solFile << x << "\t" << y << "\t" << m[i][j] << "\n";

          Real e = std::abs(m[i][j] - pot_sphere(x * x + y * y, a2 * h2));
          errFile << x << "\t" << y << "\t" << e << "\n";
          e /= std::abs(m[i][j]);
          relFile << x << "\t" << y << "\t" << e << "\n";
        }
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
  }

  delete solver;
}
