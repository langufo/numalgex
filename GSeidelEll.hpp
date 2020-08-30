#ifndef GSEIDELELL_HPP
#define GSEIDELELL_HPP

#include "EllSolver.hpp"
#include "GSeidelSolver.hpp"

class GSeidelEll
  : public GSeidelSolver
  , public EllSolver
{
public:
  GSeidelEll(const double *rhs, BndryLayout neumLayout, double x0, double y0,
             double h, int nX, int nY, const double *bottom, const double *top,
             const double *left, const double *right);
};

#endif
