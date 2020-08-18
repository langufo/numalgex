#ifndef SORPDE_HPP
#define SORPDE_HPP

#include "PDESolver.hpp"

class SORSolver : public virtual PDESolver
{
public:
  SORSolver(double (*rhs)(double, double),
            BndryLayout neumLayout,
            double x0,
            double y0,
            double h,
            int nX,
            int nY,
            const double* bottom,
            const double* top,
            const double* left,
            const double* right,
            double w);

  double iter(double* sol, bool resAsErr) final override;

private:
  double line_update(int i,
                     int jFirst,
                     int jLast,
                     BndryLayout neum,
                     const double* bottom,
                     const double* top,
                     const double* left,
                     const double* right);

  double w;
};

#endif
