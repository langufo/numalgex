#include "CylSolver.hpp"

CylSolver::CylSolver(double (*rhs)(double, double), BndryLayout neumLayout,
                     double r0, double z0, double h, int nR, int nZ,
                     const double *bottom, const double *top,
                     const double *left, const double *right)
  : PDESolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right)
{}

double
CylSolver::next_value(BndryLayout neum, double r, double z, double bottom,
                      double top, double left, double right, double rhs) const
{
  double num = 0;
  double denom = 0;

  switch (neum & (LEFTBNDRY | RIGHTBNDRY)) {
    case 0:
      denom += 6;
      num += ((right - left) * h / r / 2 + right + left) * 3;
      break;
    case LEFTBNDRY:
      denom += (1 + h / r) * 2;
      num += (right * 2 + left * h) * h / r + (right - left * h) * 2;
      break;
    case RIGHTBNDRY:
      denom += (1 - h / r) * 2;
      num += (right * h - left * 2) * h / r + (right * h + left) * 2;
      break;
    default: // LEFTBNDRY | RIGHTBNDRY
      num += ((right + left) * h / r + right - left) * h * 1.5;
  }

  switch (neum & (BOTTOMBNDRY | TOPBNDRY)) {
    case 0:
      denom += 6;
      num += (bottom + top) * 3;
      break;
    case BOTTOMBNDRY:
      denom += 2;
      num += (top - bottom * h) * 2;
      break;
    case TOPBNDRY:
      denom += 2;
      num += (bottom + top * h) * 2;
      break;
    default: // BOTTOMBNDRY | TOPBNDRY
      num += (top - bottom) * h * 1.5;
  }

  return (num - rhs * h * h * 3) / denom;
}

double
CylSolver::residual(PDESolver::BndryLayout neum, double r, double z,
                    double middle, double bottom, double top, double left,
                    double right, double rhs) const
{
  double lhs = 0;

  switch (neum & (LEFTBNDRY | RIGHTBNDRY)) {
    case 0:
      lhs += (right - left) * h / 2 / r + right + left - 2 * middle;
      break;
    case LEFTBNDRY:
      lhs += (((right - middle) * 2 + left * h) * h / r +
              (right - left * h - middle) * 2) /
             3;
      break;
    case RIGHTBNDRY:
      lhs += ((right * h + (middle - left) * 2) * h / r +
              (right * h + left - middle) * 2) /
             3;
      break;
    default: // case LEFTBNDRY | RIGHTBNDRY:
      lhs += ((right + left) * h / r + right - left) * h / 2;
  }

  switch (neum & (BOTTOMBNDRY | TOPBNDRY)) {
    case 0:
      lhs += bottom + top - middle * 2;
      break;
    case BOTTOMBNDRY:
      lhs += top + bottom * h - middle;
      break;
    case TOPBNDRY:
      lhs += bottom + top * h - middle;
      break;
    default: // case BOTTOMBNDRY | TOPBNDRY:
      lhs += (top - bottom) * h / 2;
  }

  return lhs - rhs * h * h;
}
