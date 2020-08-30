#ifndef GSEIDELPDE_HPP
#define GSEIDELPDE_HPP

#include "PDESolver.hpp"

class GSeidelSolver : public virtual PDESolver
{
public:
  GSeidelSolver(const double *rhs, BndryLayout neumLayout, double x0, double y0,
                double h, int nX, int nY, const double *bottom,
                const double *top, const double *left, const double *right);
  
  double iter(double *sol, bool resAsErr) final override;
  
  void rev_solv_direc(bool revX, bool revY);

private:
  double line_update(int i, int jInf, int jSup, BndryLayout neum,
                     const double *bottom, const double *top,
                     const double *left, const double *right);

  bool revX;
  bool revY;
};

#endif
