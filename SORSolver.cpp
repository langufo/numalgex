#include "SORSolver.hpp"

#include <math.h>

SORSolver::SORSolver(double (*rhs)(double, double),
                     PDESolver::BndryLayout neumLayout, double x0, double y0,
                     double h, int nX, int nY, const double *bottom,
                     const double *top, const double *left, const double *right,
                     double w)
  : PDESolver(rhs, neumLayout, x0, y0, h, nX, nY, bottom, top, left, right),
    w(w),
    revX(false),
    revY(false)
{}

double
SORSolver::line_update(int i, int jInf, int jSup, BndryLayout neum,
                       const double *bottom, const double *top,
                       const double *left, const double *right)
{
  double sum = 0;

  int step = revY ? -1 : 1;
  int jFirst = revY ? jSup : jInf;

  for (int j = jFirst; j >= jInf && j <= jSup; j += step) {
    double delta = w * (next_value(neum, x[i], y[j], *bottom, *top, *left,
                                   *right, mr[i][j]) -
                        ms[i][j]);

    ms[i][j] += delta;

    sum += fabs(delta);

    top += step;
    bottom += step;
    left += step;
    right += step;
  }

  return sum;
}

double
SORSolver::iter(double *sol, bool resAsErr)
{
  ms.set_first_elem(sol);

  double eps = 0;

  int hStep = revX ? -nY : nY;
  int iStep = revX ? -1 : 1;

  for (int c = 0; c < nRegX; ++c) {
    int p = revX ? nRegX - 1 - c : c;
    for (int d = 0; d < nRegY; ++d) {
      int q = revY ? nRegY - 1 - d : d;

      int iStart = revX ? limX[p][1] : limX[p][0];
      int jInf = limY[q][0];
      int jSup = limY[q][1];

      const double *b, *t, *l, *r;
      int bStep, tStep;

      if (bndryY[q] & BOTTOMBNDRY) {
        b = bottom + iStart;
        bStep = 1;
      } else {
        b = ms[iStart] + jInf - 1;
        bStep = nY;
      }
      if (bndryY[q] & TOPBNDRY) {
        t = top + iStart;
        tStep = 1;
      } else {
        t = ms[iStart] + jInf + 1;
        tStep = nY;
      }
      if (bndryX[p] & LEFTBNDRY) {
        l = left + jInf;
      } else {
        l = ms[iStart - 1] + jInf;
      }
      if (bndryX[p] & RIGHTBNDRY) {
        r = right + jInf;
      } else {
        r = ms[iStart + 1] + jInf;
      }

      if (revX) {
        bStep = -bStep;
        tStep = -tStep;
      }

      for (int i = iStart; i >= limX[p][0] && i <= limX[p][1]; i += iStep) {
        eps += line_update(i, jInf, jSup, neumLayout & (bndryX[p] | bndryY[q]),
                           b, t, l, r);
        b += bStep;
        t += tStep;
        l += hStep;
        r += hStep;
      }
    }
  }
  eps *= h * h;

  if (resAsErr) {
    eps = abs_res_sum();
  }

  return eps;
}

void
SORSolver::rev_solv_direc(bool revX, bool revY)
{
  SORSolver::revX = revX;
  SORSolver::revY = revY;
}
