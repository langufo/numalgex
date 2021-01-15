#include "JacobiSolver.hpp"

#include <cmath>
#include <vector>

#include "Matrix.hpp"
#include "Real.hpp"

JacobiSolver::JacobiSolver(Real x0, Real y0, Real h, int nX,
                           int nY)
    : PDESolver(x0, y0, h, nX, nY), bBuff(nX), lBuff(nY) {}

Real JacobiSolver::iter(Real *sol, const PDE &pde) {
  Real sum = 0; /**< errore dell'iterazione */

  /* per accedere ad array come fossero matrici */
  Matrix<Real> ms(nY, sol);
  Matrix<const Real> mr(nY, pde.rhs);

  /* inizializzo i buffer */
  for (int i = 0; i < nX; ++i) {
    bBuff[i] = pde.bottom[i];
  }
  for (int j = 0; j < nY; ++j) {
    lBuff[j] = pde.left[j];
  }

  /* scorro sulle regioni */
  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      /* limiti della regione corrente */
      int iInf = limX[p][0];
      int iSup = limX[p][1];
      int jInf = limY[q][0];
      int jSup = limY[q][1];

      const Real *t, *r; // elementi adiacenti
      int tStep;         // offset per muovere t da x a x+h
      if (bndryX[p] & RIGHTBNDRY) {
        r = pde.right + jInf;
      } else {
        r = ms[iInf + 1] + jInf;
      }
      /* l'offset da applicare a t passando da x a x+h cambia a
       * seconda che ci si trovi accanto al bordo superiore o
       * meno */
      if (bndryY[q] & TOPBNDRY) {
        t = pde.top + iInf;
        tStep = 1;
      } else {
        t = ms[iInf] + jInf + 1;
        tStep = nY;
      }

      /* determino su quali fronti è data la derivata */
      BndryProp neum = pde.neum & (bndryX[p] | bndryY[q]);

      /* scorro sui punti della regione */
      for (int i = iInf; i <= iSup; ++i) {
        for (int j = jInf; j <= jSup; ++j) {
          Real old = ms[i][j];
          int k = j - jInf;
          ms[i][j] = pde.next_value(x[i], y[j], mr[i][j], neum,
                                    bBuff[i], t[k], lBuff[j],
                                    r[k], h);
          sum += std::abs(ms[i][j] - old);

          bBuff[i] = old; // da usare poi a (i,j+1)
          lBuff[j] = old; // da usare poi a (i+1,j)
        }

        /* muovo i puntatori da x a x+h */
        t += tStep;
        r += nY;
      }
    }
  }

  return sum * h * h;
}
