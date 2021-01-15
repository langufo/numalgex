#include "CylinderPDE.hpp"

#include "BndryProp.hpp"
#include "Real.hpp"

Real CylinderPDE::next_value(Real r, Real z, Real rhs,
                             BndryProp neum, Real bottom,
                             Real top, Real left, Real right,
                             Real h) {
  Real num = 0;
  Real denom = 0;

  switch (neum & (LEFTBNDRY | RIGHTBNDRY)) {
  case 0:
    denom += r * 6;
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
  default: // case LEFTBNDRY | RIGHTBNDRY
    num += (right * (r + h) - left * (r - h)) * h * 3 / 2;
  }

  switch (neum & (BOTTOMBNDRY | TOPBNDRY)) {
  case 0:
    denom += r * 6;
    num += (bottom + top) * r * 3;
    break;
  case BOTTOMBNDRY:
    denom += r * 2;
    num += (top - bottom * h) * r * 2;
    break;
  case TOPBNDRY:
    denom += r * 2;
    num += (bottom + top * h) * r * 2;
    break;
  default: // case BOTTOMBNDRY | TOPBNDRY
    num += (top - bottom) * h * r * 3 / 2;
  }

  return (num - rhs * h * h * r * 3) / denom;
}

Real CylinderPDE::residual(Real r, Real z, Real middle,
                           Real rhs, BndryProp neum,
                           Real bottom, Real top, Real left,
                           Real right, Real h) {
  Real lhs = 0;

  switch (neum & (LEFTBNDRY | RIGHTBNDRY)) {
  case 0:
    lhs += right * (1 + h / r / 2) + left * (1 - h / r / 2) -
           middle * 2;
    break;
  case LEFTBNDRY:
    lhs += ((right - middle) * (1 + h / r) * 2 -
            left * h * (2 - h / r)) /
           3;
    break;
  case RIGHTBNDRY:
    lhs += (right * h * (2 + h / r) +
            (left - middle) * (1 - h / r) * 2) /
           3;
    break;
  default: // case LEFTBNDRY | RIGHTBNDRY:
    lhs += (right * (1 + h / r) - left * (1 - h / r)) * h / 2;
  }

  switch (neum & (BOTTOMBNDRY | TOPBNDRY)) {
  case 0:
    lhs += bottom + top - middle * 2;
    break;
  case BOTTOMBNDRY:
    lhs += (top - bottom * h - middle) * 2 / 3;
    break;
  case TOPBNDRY:
    lhs += (bottom + top * h - middle) * 2 / 3;
    break;
  default: // case BOTTOMBNDRY | TOPBNDRY:
    lhs += (top - bottom) * h / 2;
  }

  return lhs - rhs * h * h;
}
