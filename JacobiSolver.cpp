#include "JacobiSolver.hpp"

#include <math.h>
#include <vector>

#include "Matrix.hpp"

JacobiSolver::JacobiSolver(double (*rhs)(double, double),
                           PDESolver::BndryLayout l, double x0, double y0,
                           double h, int nX, int nY, const double *bottom,
                           const double *top, const double *left,
                           const double *right)
  : PDESolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right),
    hBuff(nX),
    vBuff(nY)
{}

double
JacobiSolver::line_update(int i, int jFirst, int jLast, BndryLayout neum,
                          const double *top, const double *right)
{
  double sum = 0; // sum of the absolute value of the corrections applied

  for (int j = jFirst; j <= jLast; ++j) {
    double a = ms[i][j]; // old value copy
    int k = j - jFirst;

    ms[i][j] = next_value(neum, x[i], y[j], hBuff[i], top[k], vBuff[j],
                          right[k], mr[i][j]);

    sum += fabs(ms[i][j] - a);

    hBuff[i] = a; // update the value to be used at (i, j+1)
    vBuff[j] = a; // update the value to be used at (i+1, j)
  }

  return sum;
}

double
JacobiSolver::iter(double *s, bool resAsErr)
{
  ms.set_first_elem(s);

  /* loading the buffers */
  for (int i = 0; i < nX; ++i) {
    hBuff[i] = bottom[i];
  }
  for (int j = 0; j < nY; ++j) {
    vBuff[j] = left[j];
  }

  double eps = 0;

  /*
   * two arrays are used as buffers to store the values at different points in
   * the lattice.
   * - to update the solution at (x,y), hBuff[i] has to be equal to the value at
   * (x,y-1) before it was updated, and vBuff[j] to that at (x-1,y);
   * - when the solution at (x,y) has been updated, its previous value gets
   * copied into both hBuff[i] and vBuff[j].
   * all the points in the lattice can be updated if they are visited in any
   * order that satisfies these two conditions.
   */
  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      int iFirst = limX[p][0];
      int jFirst = limY[q][0];
      int jLast = limY[q][1];

      /*
       * the pointers to the bottom and left values are not needed since they
       * are copied inside the buffers and properly updated by line_update as s
       * gets modified
       */
      const double *t, *r;
      int tStep;

      if (bndryY[q] & TOPBNDRY) {
        t = top + iFirst;
        tStep = 1;
      } else {
        t = ms[iFirst] + jFirst + 1;
        tStep = nY;
      }
      if (bndryX[p] & RIGHTBNDRY) {
        r = right + jFirst;
      } else {
        r = ms[iFirst + 1] + jFirst;
      }

      for (int i = limX[p][0]; i <= limX[p][1]; ++i) {
        eps += line_update(i, jFirst, jLast,
                           neumLayout & (bndryX[p] | bndryY[q]), t, r);
        t += tStep;
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
