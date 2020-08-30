#include "GSeidelCyl.hpp"

GSeidelCyl::GSeidelCyl(const double *rhs, PDESolver::BndryLayout neumLayout,
                       double r0, double z0, double h, int nR, int nZ,
                       const double *bottom, const double *top,
                       const double *left, const double *right)
  : PDESolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right),
    GSeidelSolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right),
    CylSolver(rhs, neumLayout, r0, z0, h, nR, nZ, bottom, top, left, right)
{}
