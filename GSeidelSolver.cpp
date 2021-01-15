#include "GSeidelSolver.hpp"

#include <cmath>
#include <vector>

#include "Matrix.hpp"
#include "Real.hpp"

GSeidelSolver::GSeidelSolver(Real x0, Real y0, Real h, int nX,
                             int nY)
    : PDESolver(x0, y0, h, nX, nY) {}

Real GSeidelSolver::iter(Real *sol, const PDE &pde) {
  Real sum = 0; /**< errore dell'iterazione */

  /* per accedere ad array come fossero matrici */
  Matrix<const Real> mr(nY, pde.rhs); // rhs
  Matrix<Real> ms(nY, sol);           // soluzione

  /* scorro sulle regioni */
  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      /* limiti della regione corrente */
      int iInf = limX[p][0];
      int iSup = limX[p][1];
      int jInf = limY[q][0];
      int jSup = limY[q][1];

      const Real *b, *t, *l, *r; // elementi adiacenti
      init_point_bndry(iInf, jInf, sol, pde, b, t, l, r);

      /* gli offset da applicare a b e t passando da x a x+h
       * cambiano a seconda che ci si trovi accanto al bordo
       * superiore o inferiore o meno */
      int bStep = bndryY[q] & BOTTOMBNDRY ? 1 : nY;
      int tStep = bndryY[q] & TOPBNDRY ? 1 : nY;

      /* determino su quali fronti è data la derivata */
      BndryProp neum = pde.neum & (bndryX[p] | bndryY[q]);

      /* scorro sui punti della regione */
      for (int i = iInf; i <= iSup; ++i) {
        for (int j = jInf; j <= jSup; ++j) {
          int k = j - jInf;
          Real old = ms[i][j];
          ms[i][j] = pde.next_value(x[i], y[j], mr[i][j], neum,
                                    b[k], t[k], l[k], r[k], h);
          sum += std::abs(ms[i][j] - old);
        }

        /* muovo i puntatori da x a x+h */
        b += bStep;
        t += tStep;
        l += nY;
        r += nY;
      }
    }
  }

  return sum * h * h;
}
