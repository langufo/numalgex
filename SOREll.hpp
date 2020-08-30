#ifndef SORELL_HPP
#define SORELL_HPP

#include "EllSolver.hpp"
#include "SORSolver.hpp"

class SOREll
  : public SORSolver
  , public EllSolver
{
public:
  SOREll(const double *rhs, BndryLayout neumLayout, double x0, double y0,
         double h, int nX, int nY, const double *bottom, const double *top,
         const double *left, const double *right, double w);
};

#endif
