#ifndef CYLSOLVER_HPP
#define CYLSOLVER_HPP

#include "PDESolver.hpp"

class CylSolver : public virtual PDESolver
{
protected:
  CylSolver(double (*rhs)(double, double), BndryLayout neumLayout, double r0,
            double z0, double h, int nR, int nZ, const double *bottom,
            const double *top, const double *left, const double *right);

  double next_value(PDESolver::BndryLayout neum, double r, double z,
                    double bottom, double top, double left, double right,
                    double rhs) const override;
  double residual(PDESolver::BndryLayout neum, double r, double z,
                  double middle, double bottom, double top, double left,
                  double right, double rhs) const override;
};

#endif
