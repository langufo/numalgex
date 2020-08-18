#ifndef GSEIDELPDE_HPP
#define GSEIDELPDE_HPP

#include "PDESolver.hpp"

class GSeidelSolver : public virtual PDESolver
{
public:
  GSeidelSolver(double (*rhs)(double, double),
                BndryLayout neumLayout,
                double x0,
                double y0,
                double h,
                int nX,
                int nY,
                const double* bottom,
                const double* top,
                const double* left,
                const double* right);
  double iter(double* sol, bool resAsErr) final;

private:
  double line_update(int i,
                     int jFirst,
                     int jLast,
                     BndryLayout neum,
                     const double* bottom,
                     const double* top,
                     const double* left,
                     const double* right);
};

#endif
