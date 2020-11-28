#ifndef CYLINDERPDE_HPP
#define CYLINDERPDE_HPP

#include "BndryProp.hpp"
#include "Real.hpp"

namespace CylinderPDE {

Real
next_value(Real r,
           Real z,
           Real rhs,
           BndryProp neum,
           Real bottom,
           Real top,
           Real left,
           Real right,
           Real h);

Real
residual(Real r,
         Real z,
         Real middle,
         Real rhs,
         BndryProp neum,
         Real bottom,
         Real top,
         Real left,
         Real right,
         Real h);

}

#endif
