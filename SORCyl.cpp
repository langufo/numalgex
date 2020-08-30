#include "SORCyl.hpp"

SORCyl::SORCyl(const double *rhs, BndryLayout neumLayout, double r0, double z0,
               double h, int nR, int nZ, const double *bottom,
               const double *top, const double *left, const double *right,
               double w)
  : PDESolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right),
    SORSolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right, w),
    CylSolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right)
{}
