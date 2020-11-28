#include "PDESolver.hpp"

#include <cmath>
#include <vector>

#include "Matrix.hpp"
#include "Real.hpp"

PDESolver::PDESolver(Real x0, Real y0, Real h, int nX, int nY)
  : h(h)
  , nX(nX)
  , nY(nY)
#ifndef ON_THE_FLY
  , x(nX)
  , y(nY)
#else
  , x(x0, h)
  , y(y0, h)
#endif
  , bndryX{ LEFTBNDRY, 0, RIGHTBNDRY }
  , bndryY{ BOTTOMBNDRY, 0, TOPBNDRY }
  , limX{ { 0, 0 }, { 1, nX - 2 } }
  , limY{ { 0, 0 }, { 1, nY - 2 } }
{

  nRegX = nX < 3 ? nX : 3;
  nRegY = nY < 3 ? nY : 3;

  /*
   * up to this point:
   * - if nRegX == 2, limX[1] is wrong;
   * - if nRegX == 3, limX[2] hasn't been initialized.
   * either way, the following piece of code solves all the aforementioned
   * problems (the same goes for limY)
   */
  limX[nRegX - 1][0] = nX - 1;
  limX[nRegX - 1][1] = nX - 1;
  limY[nRegY - 1][0] = nY - 1;
  limY[nRegY - 1][1] = nY - 1;

  /*
   * a similar problem involves bndryX and bndryY, and similarly it is solved
   */
  bndryX[nRegX - 1] |= RIGHTBNDRY;
  bndryY[nRegY - 1] |= TOPBNDRY;
}

Real
PDESolver::abs_res_sum(Real * sol, const PDE & pde) const
{
  Real sum = 0;

  Matrix<const Real> mr(nY, pde.rhs);
  Matrix<const Real> ms(nY, sol);

  /*
   * since each region has well-defined boundary properties, looping over them
   * requires less checks at runtime
   */
  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iInf = limX[p][0];
      int iSup = limX[p][1];
      int jInf = limY[q][0];
      int jSup = limY[q][1];

      BndryProp bndry = bndryX[p] | bndryY[q];

      const Real *b, *t, *l, *r; // pointers to adjacent lattice points
      init_point_bndry(iInf, jInf, sol, pde, b, t, l, r);

      /*
       * while staying in the same region, going from a point to the one on the
       * right instead requires (b,t,l,r) -> (b+nY,t+nY,l+nY,r+nY), UNLESS b (t)
       * is on the bottom or top boundary: then b -> b+1 or t -> t+1
       */
      int bStep = bndry & BOTTOMBNDRY ? 1 : nY; // x -> x+h
      int tStep = bndry & TOPBNDRY ? 1 : nY;    // x -> x+h

      for (int i = iInf; i <= iSup; ++i) {
        for (int j = jInf; j <= jSup; ++j) {
          int k = j - jInf;
          sum += std::fabs(pde.residual(x[i],
                                        y[j],
                                        ms[i][j],
                                        mr[i][j],
                                        pde.neum & bndry,
                                        b[k],
                                        t[k],
                                        l[k],
                                        r[k],
                                        h));
        }

        /* moving the pointers from x to x+h */
        b += bStep;
        t += tStep;
        l += nY;
        r += nY;
      }
    }
  }

  return sum;
}

void
PDESolver::init_point_bndry(int i,
                            int j,
                            const Real * sol,
                            const PDE & pde,
                            const Real *& b,
                            const Real *& t,
                            const Real *& l,
                            const Real *& r) const
{
  if (j == 0) {
    b = pde.bottom + i;
  } else {
    b = sol + i * nY + j - 1;
  }

  if (j == nY - 1) {
    t = pde.top + i;
  } else {
    t = sol + i * nY + j + 1;
  }

  if (i == 0) {
    l = pde.left + j;
  } else {
    l = sol + (i - 1) * nY + j;
  }

  if (i == nX - 1) {
    r = pde.right + j;
  } else {
    r = sol + (i + 1) * nY + j;
  }
}
