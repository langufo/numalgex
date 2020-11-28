#ifndef GSEIDELPDE_HPP
#define GSEIDELPDE_HPP

#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"

class GSeidelSolver : public PDESolver
{
public:
  GSeidelSolver(Real x0, Real y0, Real h, int nX, int nY);

  Real iter(Real * sol, const PDE & pde) override;
};

#endif
