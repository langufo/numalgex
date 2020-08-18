#include "JacobiCyl.hpp"

JacobiCyl::JacobiCyl(double (*rhs)(double, double),
                     PDESolver::BndryLayout l,
                     double r0,
                     double z0,
                     double h,
                     int nR,
                     int nZ,
                     const double* bottom,
                     const double* top,
                     const double* left,
                     const double* right)
  : PDESolver(rhs, l, r0, z0, h, nR, nZ, bottom, top, left, right)
  , JacobiSolver(rhs, l, r0, z0, h, nR, nZ, bottom, top, left, right)
  , CylSolver(rhs, l, r0, z0, h, nR, nZ, bottom, top, left, right)
{}
