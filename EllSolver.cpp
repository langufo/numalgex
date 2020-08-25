#include "EllSolver.hpp"

EllSolver::EllSolver(double (*rhs)(double, double), BndryLayout neumLayout,
                     double x0, double y0, double h, int nX, int nY,
                     const double *bottom, const double *top,
                     const double *left, const double *right)
  : PDESolver(rhs, neumLayout, x0, y0, h, nX, nY, bottom, top, left, right)
{}

double
EllSolver::next_value(BndryLayout neum, double x, double y, double bottom,
                      double top, double left, double right, double rhs) const
{
  double a = 0;
  double n = 0;

  BndryLayout mask[2][2] = { { BOTTOMBNDRY, TOPBNDRY },
                             { LEFTBNDRY, RIGHTBNDRY } };
  double val[2][2] = { { bottom, top }, { left, right } };

  for (int k = 0; k < 2; ++k) {
    if (bool(neum & mask[k][0]) == bool(neum & mask[k][1])) {
      if (neum & mask[k][0]) {
        a += (val[k][1] - val[k][0]) * h / 2;
      } else {
        n += 2;
        a += val[k][0] + val[k][1];
      }
    } else {
      n += 2.0 / 3.0;
      if (neum & mask[k][0]) {
        a += (val[k][1] - val[k][0] * h) * 2 / 3;
      } else {
        a += (val[k][0] + val[k][1] * h) * 2 / 3;
      }
    }
  }

  return (a - rhs * h * h) / n;
}

double
EllSolver::residual(BndryLayout neum, double x, double y, double middle,
                    double bottom, double top, double left, double right,
                    double rhs) const
{
  double a = 0;

  BndryLayout mask[2][2] = { { BOTTOMBNDRY, TOPBNDRY },
                             { LEFTBNDRY, RIGHTBNDRY } };
  double val[2][2] = { { bottom, top }, { left, right } };

  for (int k = 0; k < 2; ++k) {
    if (bool(neum & mask[k][0]) == bool(neum & mask[k][1])) {
      if (neum & mask[k][0]) {
        a += (val[k][1] - val[k][0]) * h / 2;
      } else {
        a += val[k][0] + val[k][1] - middle * 2;
      }
    } else {
      if (neum & mask[k][0]) {
        a += (val[k][1] - val[k][0] * h - middle) * 2 / 3;
      } else {
        a += (val[k][0] + val[k][1] * h - middle) * 2 / 3;
      }
    }
  }

  return a - rhs * h * h;
}
