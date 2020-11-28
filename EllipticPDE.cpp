#include "EllipticPDE.hpp"

#include "BndryProp.hpp"
#include "Real.hpp"

Real
EllipticPDE::next_value(Real x,
                        Real y,
                        Real rhs,
                        Real h,
                        BndryProp neum,
                        Real bottom,
                        Real top,
                        Real left,
                        Real right)
{
  Real a = 0;
  Real n = 0;

  BndryProp mask[2][2] = { { BOTTOMBNDRY, TOPBNDRY },
                           { LEFTBNDRY, RIGHTBNDRY } };
  Real val[2][2] = { { bottom, top }, { left, right } };

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

Real
EllipticPDE::residual(Real x,
                      Real y,
                      Real middle,
                      Real rhs,
                      Real h,
                      BndryProp neum,
                      Real bottom,
                      Real top,
                      Real left,
                      Real right)
{
  Real a = 0;

  BndryProp mask[2][2] = { { BOTTOMBNDRY, TOPBNDRY },
                           { LEFTBNDRY, RIGHTBNDRY } };
  Real val[2][2] = { { bottom, top }, { left, right } };

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
