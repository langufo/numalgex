#include "PDESolver.hpp"

#include <math.h>

#include <vector>

#include "Matrix.hpp"

PDESolver::PDESolver(double (*rhs)(double, double),
                     BndryLayout neumLayout,
                     double x0,
                     double y0,
                     double h,
                     int nX,
                     int nY,
                     const double* bottom,
                     const double* top,
                     const double* left,
                     const double* right)
  : rhs(rhs)
  , h(h)
  , nX(nX)
  , nY(nY)
  , x(nX)
  , y(nY)
  , neumLayout(neumLayout)
  , bottom(bottom)
  , top(top)
  , left(left)
  , right(right)
  , r(nX * nY)
  , bndryX{ NOBNDRY, NOBNDRY, NOBNDRY }
  , bndryY{ NOBNDRY, NOBNDRY, NOBNDRY }
  , limX{ { 0, 0 }, { 1, nX - 2 } }
  , limY{ { 0, 0 }, { 1, nY - 2 } }
  , mr(nY, r.data())
  , ms(nY)
{
  for (int i = 0; i < nX; ++i) {
    x[i] = x0 + i * h;
  }
  for (int j = 0; j < nY; ++j) {
    y[j] = y0 + j * h;
  }

  nRegX = nX < 3 ? nX : 3;
  nRegY = nY < 3 ? nY : 3;

  limX[nRegX - 1][0] = nX - 1;
  limX[nRegX - 1][1] = nX - 1;
  limY[nRegY - 1][0] = nY - 1;
  limY[nRegY - 1][1] = nY - 1;

  bndryX[0] |= LEFTBNDRY;
  bndryX[nRegX - 1] |= RIGHTBNDRY;

  bndryY[0] |= BOTTOMBNDRY;
  bndryY[nRegY - 1] |= TOPBNDRY;

  cache_rhs();
}

void
PDESolver::cache_rhs()
{
  for (int i = 0; i < nX; ++i) {
    for (int j = 0; j < nY; ++j) {
      mr[i][j] = rhs(x[i], y[j]);
    }
  }
}

double
PDESolver::line_error(int i,
                      int jFirst,
                      int jLast,
                      BndryLayout neum,
                      const double* bottom,
                      const double* top,
                      const double* left,
                      const double* right)
{
  double sum = 0;

  for (int j = jFirst; j <= jLast; ++j) {
    int k = j - jFirst;
    sum += fabs(residual(neum,
                         x[i],
                         y[j],
                         ms[i][j],
                         bottom[k],
                         top[k],
                         left[k],
                         right[k],
                         mr[i][j]));
  }

  return sum;
}

double
PDESolver::abs_res_sum()
{
  double sum = 0;

  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iFirst = limX[p][0];
      int jFirst = limY[q][0];
      int jLast = limY[q][1];

      BndryLayout bndry = bndryX[p] | bndryY[q];

      const double *b, *t, *l, *r;
      int bStep, tStep;

      if (bndry & BOTTOMBNDRY) {
        b = bottom + iFirst;
        bStep = 1;
      } else {
        b = ms[iFirst] + jFirst - 1;
        bStep = nY;
      }
      if (bndry & TOPBNDRY) {
        t = top + iFirst;
        tStep = 1;
      } else {
        t = ms[iFirst] + jFirst + 1;
        tStep = nY;
      }
      if (bndry & LEFTBNDRY) {
        l = left + jFirst;
      } else {
        l = ms[iFirst - 1] + jFirst;
      }
      if (bndry & RIGHTBNDRY) {
        r = right + jFirst;
      } else {
        r = ms[iFirst + 1] + jFirst;
      }

      for (int i = iFirst; i <= limX[p][1]; ++i) {
        sum += line_error(i, jFirst, jLast, neumLayout & bndry, b, t, l, r);

        b += bStep;
        t += tStep;
        l += nY;
        r += nY;
      }
    }
  }

  return sum;
}
