#ifndef SORCYL_HPP
#define SORCYL_HPP

#include "CylSolver.hpp"
#include "SORSolver.hpp"

class SORCyl
  : public SORSolver
  , public CylSolver
{
public:
  SORCyl(double (*rhs)(double, double),
         BndryLayout neumLayout,
         double r0,
         double z0,
         double h,
         int nR,
         int nZ,
         const double* bottom,
         const double* top,
         const double* left,
         const double* right,
         double w);
};

#endif
