#ifndef ELLSOLVER_HPP
#define ELLSOLVER_HPP

#include "PDESolver.hpp"

class EllSolver : public virtual PDESolver
{
protected:
  EllSolver(double (*rhs)(double, double), BndryLayout neumLayout, double x0,
            double y0, double h, int nX, int nY, const double *bottom,
            const double *top, const double *left, const double *right);

  double next_value(PDESolver::BndryLayout neum, double x, double y,
                    double bottom, double top, double left, double right,
                    double rhs) const override;
  double residual(PDESolver::BndryLayout neum, double x, double y,
                  double middle, double bottom, double top, double left,
                  double right, double rhs) const override;
};

#endif
