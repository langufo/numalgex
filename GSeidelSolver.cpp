#include "GSeidelSolver.hpp"

#include <cmath>
#include <vector>

#include "Matrix.hpp"
#include "Real.hpp"

GSeidelSolver::GSeidelSolver(Real x0, Real y0, Real h, int nX, int nY)
  : PDESolver(x0, y0, h, nX, nY)
{}

Real
GSeidelSolver::iter(Real * sol, const PDE & pde)
{
  Real sum = 0;

  Matrix<const Real> mr(nY, pde.rhs);
  Matrix<Real> ms(nY, sol);

  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iInf = limX[p][0];
      int iSup = limX[p][1];
      int jInf = limY[q][0];
      int jSup = limY[q][1];

      const Real *b, *t, *l, *r;

      init_point_bndry(iInf, jInf, sol, pde, b, t, l, r);

      int bStep = bndryY[q] & BOTTOMBNDRY ? 1 : nY;
      int tStep = bndryY[q] & TOPBNDRY ? 1 : nY;

      BndryProp neum = pde.neum & (bndryX[p] | bndryY[q]);

      for (int i = iInf; i <= iSup; ++i) {
        for (int j = jInf; j <= jSup; ++j) {
          Real old = ms[i][j];

          int k = j - jInf;

          ms[i][j] = pde.next_value(
            x[i], y[j], mr[i][j], neum, b[k], t[k], l[k], r[k], h);

          sum += std::fabs(ms[i][j] - old);
        }

        b += bStep;
        t += tStep;
        l += nY;
        r += nY;
      }
    }
  }

  return sum * h * h;
}
