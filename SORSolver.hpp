#ifndef SORPDE_HPP
#define SORPDE_HPP

#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"

class SORSolver : public PDESolver
{
public:
  SORSolver(Real x0, Real y0, Real h, int nX, int nY, Real w);

  Real iter(Real * sol, const PDE & pde) override;

  void rev_solv_direc(bool revX, bool revY);

private:
  Real w;
};

#endif
