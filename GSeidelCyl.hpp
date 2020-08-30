#ifndef GSEIDELCYL_HPP
#define GSEIDELCYL_HPP

#include "CylSolver.hpp"
#include "GSeidelSolver.hpp"

class GSeidelCyl
  : public GSeidelSolver
  , public CylSolver
{
public:
  GSeidelCyl(const double *rhs, BndryLayout neumLayout, double r0, double z0,
             double h, int nR, int nZ, const double *bottom, const double *top,
             const double *left, const double *right);
};

#endif
