#include "CylSolver.hpp"

CylSolver::CylSolver(double (*rhs)(double, double),
                     PDESolver::BndryLayout neumLayout,
                     double r0,
                     double z0,
                     double h,
                     int nR,
                     int nZ,
                     const double* bottom,
                     const double* top,
                     const double* left,
                     const double* right)
  : PDESolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right)
{}

double
CylSolver::next_value(PDESolver::BndryLayout neum,
                      double r,
                      double z,
                      double bottom,
                      double top,
                      double left,
                      double right,
                      double rhs)
{
  double a = 0;
  double c = 0;

  if (bool(neum & LEFTBNDRY) == bool(neum & RIGHTBNDRY)) {
    if (neum & LEFTBNDRY) {
      a += ((right + left) * h / r + right - left) * h / 2;
    } else {
      c += 2;
      a += (right - left) * h / 2 / r + right + left;
    }
  } else {
    if (neum & LEFTBNDRY) {
      c += 2.0 / 3.0 * (1 + h / r);
      a += ((right * 2 + left * h) * h / r + (right - left * h) * 2) / 3;
    } else {
      c += 2.0 / 3.0 * (1 - h / r);
      a += ((right * h - left * 2) * h / r + (right * h + left) * 2) / 3;
    }
  }

  if (bool(neum & BOTTOMBNDRY) == bool(neum & TOPBNDRY)) {
    if (neum & BOTTOMBNDRY) {
      a += (top - bottom) * h / 2;
    } else {
      c += 2;
      a += bottom + top;
    }
  } else {
    c += 2.0 / 3.0;
    if (neum & BOTTOMBNDRY) {
      a += (top - bottom * h) * 2 / 3;
    } else {
      a += (bottom + top * h) * 2 / 3;
    }
  }

  return (a - rhs * h * h) / c;
}

double
CylSolver::residual(PDESolver::BndryLayout neum,
                    double r,
                    double z,
                    double middle,
                    double bottom,
                    double top,
                    double left,
                    double right,
                    double rhs)
{
  double a = 0;

  if (bool(neum & LEFTBNDRY) == bool(neum & RIGHTBNDRY)) {
    if (neum & LEFTBNDRY) {
      a += ((right + left) * h / r + right - left) * h / 2;
    } else {
      a += (right - left) * h / 2 / r + right + left - 2 * middle;
    }
  } else {
    if (neum & LEFTBNDRY) {
      a += (((right - middle) * 2 + left * h) * h / r +
            (right - left * h - middle) * 2) /
           3;
    } else {
      a += ((right * h + (middle - left) * 2) * h / r +
            (right * h + left - middle) * 2) /
           3;
    }
  }

  if (bool(neum & BOTTOMBNDRY) == bool(neum & TOPBNDRY)) {
    if (neum & BOTTOMBNDRY) {
      a += (top - bottom) * h / 2;
    } else {
      a += bottom + top - middle * 2;
    }
  } else {
    if (neum & BOTTOMBNDRY) {
      a += top + bottom * h - middle;
    } else {
      a += bottom + top * h - middle;
    }
  }

  return a - rhs * h * h;
}
