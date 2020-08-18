#ifndef JACOBICYL_HPP
#define JACOBICYL_HPP

#include "CylSolver.hpp"
#include "JacobiSolver.hpp"

class JacobiCyl
  : public JacobiSolver
  , public CylSolver
{
public:
  JacobiCyl(double (*rhs)(double, double),
            BndryLayout l,
            double r0,
            double z0,
            double h,
            int nR,
            int nZ,
            const double* bottom,
            const double* top,
            const double* left,
            const double* right);
};

#endif
