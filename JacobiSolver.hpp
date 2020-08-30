#ifndef JACOBIPDE_HPP
#define JACOBIPDE_HPP

#include <vector>

#include "PDESolver.hpp"

class JacobiSolver : public virtual PDESolver
{
public:
  JacobiSolver(const double *rhs, BndryLayout l, double x0, double y0, double h,
               int nX, int nY, const double *bottom, const double *top,
               const double *left, const double *right);

  double iter(double *s, bool resAsErr) final;

private:
  double line_update(int i, int jFirst, int jLast, BndryLayout neum,
                     const double *top, const double *right);

  std::vector<double> hBuff, vBuff;
};

#endif
