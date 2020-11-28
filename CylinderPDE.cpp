#include "CylinderPDE.hpp"

#include "BndryProp.hpp"
#include "Real.hpp"

#define BETTER

Real
CylinderPDE::next_value(Real r,
                        Real z,
                        Real rhs,
                        BndryProp neum,
                        Real bottom,
                        Real top,
                        Real left,
                        Real right,
                        Real h)
{
  Real num = 0;
  Real denom = 0;

  switch (neum & (LEFTBNDRY | RIGHTBNDRY)) {
    case 0:
      denom += 6 * r;
      num += (right * (r + h / 2) + left * (r - h / 2)) * 3;
      break;
    case LEFTBNDRY:
      denom += (r + h) * 2;
      num += right * (r + h) * 2 - left * h * (2 * r - h);
      break;
    case RIGHTBNDRY:
      denom += (r - h) * 2;
      num += right * h * (2 * r + h) + left * (r - h) * 2;
      break;
    default: // LEFTBNDRY | RIGHTBNDRY
      num += (right * (r + h) - left * (r - h)) * h * (static_cast<Real>(3)/2);
  }

  switch (neum & (BOTTOMBNDRY | TOPBNDRY)) {
    case 0:
      denom += 6 * r;
      num += (bottom + top) * 3 * r;
      break;
    case BOTTOMBNDRY:
      denom += 2 * r;
      num += (top - bottom * h) * 2 * r;
      break;
    case TOPBNDRY:
      denom += 2 * r;
      num += (bottom + top * h) * 2 * r;
      break;
    default: // BOTTOMBNDRY | TOPBNDRY
      num += (top - bottom) * h * (static_cast<Real>(3)/2) * r;
  }

  return (num - rhs * h * h * 3 * r) / denom;
}

Real
CylinderPDE::residual(Real r,
                      Real z,
                      Real middle,
                      Real rhs,
                      BndryProp neum,
                      Real bottom,
                      Real top,
                      Real left,
                      Real right,
                      Real h)
{
  Real lhs = 0;

#ifndef BETTER
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
#else
  switch (neum & (LEFTBNDRY | RIGHTBNDRY)) {
    case 0:
      lhs += right * (1 + h / r / 2) + left * (1 - h / r / 2) - 2 * middle;
      break;
    case LEFTBNDRY:
      lhs += ((right - middle) * (1 + h / r) * 2 - left * h * (2 - h / r)) / 3;
      break;
    case RIGHTBNDRY:
      lhs += (right * h * (2 + h / r) + (left - middle) * (1 - h / r) * 2) / 3;
      break;
    default: // case LEFTBNDRY | RIGHTBNDRY:
      lhs += (right * (1 + h / r) - left * (1 - h / r)) * h / 2;
  }
#endif

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
