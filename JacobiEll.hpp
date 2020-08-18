#ifndef JACOBIELL_HPP
#define JACOBIELL_HPP

#include "EllSolver.hpp"
#include "JacobiSolver.hpp"

class JacobiEll
  : public JacobiSolver
  , public EllSolver
{
public:
  JacobiEll(double (*rhs)(double, double),
            BndryLayout l,
            double x0,
            double y0,
            double h,
            int nX,
            int nY,
            const double* bottom,
            const double* top,
            const double* left,
            const double* right);
};

#endif
