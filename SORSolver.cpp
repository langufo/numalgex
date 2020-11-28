#include "SORSolver.hpp"

#include <cmath>

#include "Matrix.hpp"
#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"

SORSolver::SORSolver(Real x0, Real y0, Real h, int nX, int nY, Real w)
  : PDESolver(x0, y0, h, nX, nY)
  , w(w)
{}

Real
SORSolver::iter(Real * sol, const PDE & pde)
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
          int k = j - jInf;
          Real a = pde.next_value(
            x[i], y[j], mr[i][j], neum, b[k], t[k], l[k], r[k], h);

          Real delta = w * (a - ms[i][j]);
          sum += std::fabs(delta);

          ms[i][j] = (1 - w) * ms[i][j] + w * a;
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
