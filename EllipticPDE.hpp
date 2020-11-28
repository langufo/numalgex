#ifndef ELLIPTICPDE_HPP
#define ELLIPTICPDE_HPP

#include "BndryProp.hpp"
#include "Real.hpp"

namespace EllipticPDE {

Real
next_value(Real r,
           Real z,
           Real rhs,
           Real h,
           BndryProp neum,
           Real bottom,
           Real top,
           Real left,
           Real right);

Real
residual(Real r,
         Real z,
         Real middle,
         Real rhs,
         Real h,
         BndryProp neum,
         Real bottom,
         Real top,
         Real left,
         Real right);

}

#endif
