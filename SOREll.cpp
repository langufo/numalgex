#include "SOREll.hpp"

SOREll::SOREll(const double *rhs, BndryLayout neumLayout, double x0, double y0,
               double h, int nX, int nY, const double *bottom,
               const double *top, const double *left, const double *right,
               double w)
  : PDESolver(rhs, neumLayout, x0, y0, h, nX, nY, bottom, top, left, right),
    SORSolver(rhs, neumLayout, x0, y0, h, nX, nY, bottom, top, left, right, w),
    EllSolver(rhs, neumLayout, x0, y0, h, nX, nY, bottom, top, left, right)
{}
