#include "SORSolver.hpp"

#include <math.h>

SORSolver::SORSolver(double (*rhs)(double, double),
                     PDESolver::BndryLayout neumLayout,
                     double x0,
                     double y0,
                     double h,
                     int nX,
                     int nY,
                     const double* bottom,
                     const double* top,
                     const double* left,
                     const double* right,
                     double w)
  : PDESolver(rhs, neumLayout, x0, y0, h, nX, nY, bottom, top, left, right)
  , w(w)
{}

double
SORSolver::line_update(int i,
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
    double a =
      w * (next_value(
             neum, x[i], y[j], bottom[k], top[k], left[k], right[k], mr[i][j]) -
           ms[i][j]);
    sum += fabs(a);
    ms[i][j] += a;
  }

  return sum;
}

double
SORSolver::iter(double* sol, bool resAsErr)
{
  ms.set_first_elem(sol);

  double eps = 0;

  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iFirst = limX[p][0];
      int jFirst = limY[q][0];
      int jLast = limY[q][1];

      const double *b, *t, *l, *r;
      int bStep, tStep;

      if (bndryY[q] & BOTTOMBNDRY) {
        b = bottom + iFirst;
        bStep = 1;
      } else {
        b = ms[iFirst] + jFirst - 1;
        bStep = nY;
      }
      if (bndryY[q] & TOPBNDRY) {
        t = top + iFirst;
        tStep = 1;
      } else {
        t = ms[iFirst] + jFirst + 1;
        tStep = nY;
      }
      if (bndryX[p] & LEFTBNDRY) {
        l = left + jFirst;
      } else {
        l = ms[iFirst - 1] + jFirst;
      }
      if (bndryX[p] & RIGHTBNDRY) {
        r = right + jFirst;
      } else {
        r = ms[iFirst + 1] + jFirst;
      }

      for (int i = limX[p][0]; i <= limX[p][1]; ++i) {
        eps += line_update(
          i, jFirst, jLast, neumLayout & (bndryX[p] | bndryY[q]), b, t, l, r);
        b += bStep;
        t += tStep;
        l += nY;
        r += nY;
      }
    }
  }
  eps *= h * h;

  if (resAsErr) {
    eps = abs_res_sum();
  }

  return eps;
}
