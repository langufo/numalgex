#ifndef JACOBISOLVER_HPP
#define JACOBISOLVER_HPP

#include <vector>

#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"

class JacobiSolver : public PDESolver {
public:
  JacobiSolver(Real x0, Real y0, Real h, int nX, int nY);

  Real iter(Real *sol, const PDE &pde) override;

private:
  /**
   * Vettore che ospita valori della soluzione precedente
   */
  std::vector<Real> bBuff;

  /**
   * Vettore che ospita valori della soluzione precedente
   */
  std::vector<Real> lBuff;
};

#endif
