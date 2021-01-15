#ifndef SORSOLVER_HPP
#define SORSOLVER_HPP

#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"

class SORSolver : public PDESolver {
public:
  SORSolver(Real x0, Real y0, Real h, int nX, int nY, Real w);

  Real iter(Real *sol, const PDE &pde) override;

private:
  Real w; /**< Parametro di sovrarilassamento */
};

#endif
