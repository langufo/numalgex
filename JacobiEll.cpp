#include "JacobiEll.hpp"

JacobiEll::JacobiEll(const double *rhs, PDESolver::BndryLayout l, double x0,
                     double y0, double h, int nX, int nY, const double *bottom,
                     const double *top, const double *left, const double *right)
  : PDESolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right),
    JacobiSolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right),
    EllSolver(rhs, l, x0, y0, h, nX, nY, bottom, top, left, right)
{}
