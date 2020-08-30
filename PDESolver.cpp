#include "PDESolver.hpp"

#include <math.h>

#include <vector>

#include "Matrix.hpp"

PDESolver::PDESolver(const double *rhs, BndryLayout neumLayout, double x0,
                     double y0, double h, int nX, int nY, const double *bottom,
                     const double *top, const double *left, const double *right)
  : h(h),
    nX(nX),
    nY(nY),
    x(nX),
    y(nY),
    neumLayout(neumLayout),
    bottom(bottom, bottom + nX),
    top(top, top + nX),
    left(left, left + nY),
    right(right, right + nY),
    bndryX{ LEFTBNDRY, 0, RIGHTBNDRY },
    bndryY{ BOTTOMBNDRY, 0, TOPBNDRY },
    limX{ { 0, 0 }, { 1, nX - 2 } },
    limY{ { 0, 0 }, { 1, nY - 2 } },
    ms(nY),
    r(rhs, rhs + nX * nY),
    mr(nY, r.data())
{
  for (int i = 0; i < nX; ++i) {
    x[i] = x0 + i * h;
  }
  for (int j = 0; j < nY; ++j) {
    y[j] = y0 + j * h;
  }

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

double
PDESolver::line_error(int i, int jFirst, int jLast, BndryLayout neum,
                      const double *bottom, const double *top,
                      const double *left, const double *right) const
{
  double sum = 0;

  for (int j = jFirst; j <= jLast; ++j) {
    int k = j - jFirst;
    sum += fabs(residual(neum, x[i], y[j], ms[i][j], bottom[k], top[k], left[k],
                         right[k], mr[i][j]));
  }

  return sum;
}

double
PDESolver::abs_res_sum() const
{
  double sum = 0;

  /*
   * since each region has well-defined boundary properties, looping over them
   * avoids a lot of checks at runtime
   */
  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iFirst = limX[p][0];
      int jFirst = limY[q][0];
      int jLast = limY[q][1];

      BndryLayout bndry = bndryX[p] | bndryY[q];

      const double *b, *t, *l, *r; // pointers to adjacent lattice points
      int bStep, tStep; // offsets to apply going from (x,y) to (x+h,y)

      /*
       *      t0   t1   t2
       *    +----+----+----+
       * l2 | s2 | s5 | s8 | r2
       *    +----+----+----+
       * l1 | s1 | s4 | s7 | r1
       *    +----+----+----+
       * l0 | s0 | s3 | s6 | r0
       *    +----+----+----+
       *      b0   b1   b2
       * going from a point to the one above (movement along the y axis)
       * requires (b,t,l,r) -> (b+1,t+1,l+1,r+1); going to the one on the right
       * instead requires (b,t,l,r) -> (b+nY,t+nY,l+nY,r+nY), UNLESS b (t) is on
       * the bottom (top) boundary: then b -> b+1 (t -> t+1)
       */
      if (bndry & BOTTOMBNDRY) {
        b = bottom.data() + iFirst;
        bStep = 1;
      } else {
        b = ms[iFirst] + jFirst - 1;
        bStep = nY;
      }
      if (bndry & TOPBNDRY) {
        t = top.data() + iFirst;
        tStep = 1;
      } else {
        t = ms[iFirst] + jFirst + 1;
        tStep = nY;
      }
      if (bndry & LEFTBNDRY) {
        l = left.data() + jFirst;
      } else {
        l = ms[iFirst - 1] + jFirst;
      }
      if (bndry & RIGHTBNDRY) {
        r = right.data() + jFirst;
      } else {
        r = ms[iFirst + 1] + jFirst;
      }

      for (int i = iFirst; i <= limX[p][1]; ++i) {
        sum += line_error(i, jFirst, jLast, neumLayout & bndry, b, t, l, r);

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
