#include "JacobiSolver.hpp"

#include <cmath>
#include <vector>

#include "Matrix.hpp"
#include "Real.hpp"

JacobiSolver::JacobiSolver(Real x0, Real y0, Real h, int nX, int nY)
  : PDESolver(x0, y0, h, nX, nY)
  , hBuff(nX)
  , vBuff(nY)
{}

Real
JacobiSolver::iter(Real * sol, const PDE & pde)
{
  Real sum = 0;
  Matrix<Real> ms(nY, sol);
  Matrix<const Real> mr(nY, pde.rhs);

  /*
   * two arrays are used as buffers to store the values at different points in
   * the lattice.
   * - to update the solution at (x,y), hBuff[i] has to be equal to the value at
   * (x,y-1) before it was updated, and vBuff[j] to the old value at (x-1,y);
   * - when the solution at (x,y) has been updated, its previous value gets
   * copied into both hBuff[i] and vBuff[j].
   * all the points in the lattice can be updated if they are visited in any
   * order that satisfies these two conditions.
   */

  /* loading the buffers */
  for (int i = 0; i < nX; ++i) {
    hBuff[i] = pde.bottom[i];
  }
  for (int j = 0; j < nY; ++j) {
    vBuff[j] = pde.left[j];
  }

  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iFirst = limX[p][0];
      int jFirst = limY[q][0];
      int jLast = limY[q][1];

      /*
       * the pointers to the bottom and left values are not needed since they
       * are copied inside the buffers and properly updated as the solution
       * changes
       */
      const Real *t, *r;
      int tStep;

      if (bndryY[q] & TOPBNDRY) {
        t = pde.top + iFirst;
        tStep = 1;
      } else {
        t = ms[iFirst] + jFirst + 1;
        tStep = nY;
      }
      if (bndryX[p] & RIGHTBNDRY) {
        r = pde.right + jFirst;
      } else {
        r = ms[iFirst + 1] + jFirst;
      }

      BndryProp neum = pde.neum & (bndryX[p] | bndryY[q]);

      for (int i = limX[p][0]; i <= limX[p][1]; ++i) {
        for (int j = jFirst; j <= jLast; ++j) {
          Real old = ms[i][j]; // old value copy
          int k = j - jFirst;

          ms[i][j] = pde.next_value(
            x[i], y[j], mr[i][j], neum, hBuff[i], t[k], vBuff[j], r[k], h);

          sum += std::fabs(ms[i][j] - old);

          hBuff[i] = old; // set up for use at (i, j+1)
          vBuff[j] = old; // set up for use at (i+1, j)
        }

        t += tStep;
        r += nY;
      }
    }
  }

  return sum * h * h;
}
