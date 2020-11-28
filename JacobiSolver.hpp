#ifndef JACOBIPDE_HPP
#define JACOBIPDE_HPP

#include <vector>

#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"

class JacobiSolver : public PDESolver
{
public:
  JacobiSolver(Real x0, Real y0, Real h, int nX, int nY);

  Real iter(Real * sol, const PDE & pde) override;

private:
  /**
   * Vector spanning along the x axis.
   */
  std::vector<Real> hBuff;

  /**
   * Vector spanning along the y axis.
   */
  std::vector<Real> vBuff;
};

#endif
