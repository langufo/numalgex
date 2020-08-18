#include "GSeidelEll.hpp"

GSeidelEll::GSeidelEll(double (*rhs)(double, double),
                       PDESolver::BndryLayout l,
                       double x0,
                       double y0,
                       double h,
                       int nX,
                       int nY,
                       const double* bottom,
                       const double* top,
                       const double* left,
                       const double* right)
  : PDESolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right)
  , GSeidelSolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right)
  , EllSolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right)
{}
