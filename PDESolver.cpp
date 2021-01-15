#include "PDESolver.hpp"

#include <cmath>
#include <vector>

#include "Real.hpp"

PDESolver::PDESolver(Real x0, Real y0, Real h, int nX, int nY)
    : h(h), nX(nX), nY(nY), x(nX),
      y(nY), bndryX{LEFTBNDRY, 0, RIGHTBNDRY},
      bndryY{BOTTOMBNDRY, 0, TOPBNDRY},
      limX{{0, 0}, {1, nX - 2}}, limY{{0, 0}, {1, nY - 2}} {
  for (int i = 0; i < nX; ++i) {
    x[i] = x0 + i * h;
  }
  for (int j = 0; j < nY; ++j) {
    y[j] = y0 + j * h;
  }

  /*
   * per la direzione x:
   * - se n == 1 c'è solo una regione, la quale si affaccia sia
   * sul bordo sinistro che su quello destro;
   * - se n == 2 ci sono due regioni, una si affaccia solo sul
   * bordo sinistro, l'altra solo su quello destro
   * - se n >= 3 ci sono tre regioni, una affacciata sul bordo
   * sinistro, una affacciata sul bordo destro e una centrale
   * che non si affaccia su nessuno dei due analogo discorso
   * vale per la direzione y
   */
  nRegX = nX < 3 ? nX : 3;
  nRegY = nY < 3 ? nY : 3;

  /*
   * fino a questo punto:
   * - se nReg* == 2, lim*[1] contiene il valore sbagliato;
   * - se nReg* == 3, lim*[2] non è stato inizializzato;
   * in ogni caso, il codice seguente risolve il problema
   */
  limX[nRegX - 1][0] = nX - 1;
  limX[nRegX - 1][1] = nX - 1;
  limY[nRegY - 1][0] = nY - 1;
  limY[nRegY - 1][1] = nY - 1;

  /*
   * un problema analogo al precedente affligge bndryX e bndrY,
   * e analogamente viene risolto
   */
  bndryX[nRegX - 1] |= RIGHTBNDRY;
  bndryY[nRegY - 1] |= TOPBNDRY;
}

Real PDESolver::abs_res_sum(Real *sol, const PDE &pde) const {
  Real sum = 0; // somma residuale

  /*
   * dal momento che le regioni hanno bordi con proprietà ben
   * definite, scorrere su di loro richiede meno controlli
   * durante l'esecuzione
   */
  for (int p = 0; p < nRegX; ++p) {
    for (int q = 0; q < nRegY; ++q) {
      /* limiti della regione corrente */
      int iInf = limX[p][0];
      int iSup = limX[p][1];
      int jInf = limY[q][0];
      int jSup = limY[q][1];

      /* bordi adiacenti alla regione */
      BndryProp bndry = bndryX[p] | bndryY[q];

      const Real *b, *t, *l, *r; // elementi adiacenti
      init_point_bndry(iInf, jInf, sol, pde, b, t, l, r);

      /*
       * restando all'interno di una regione, spostarsi da x a
       * x+h richiede (b,t,l,r)->(b+nY,t+nY,l+nY,r+nY), TRANNE
       * quando b (o t) punta al bordo inferiore (o superiore):
       * in tal caso b->b+1 (o t->t+1)
       */
      int bStep = bndry & BOTTOMBNDRY ? 1 : nY; // x -> x+h
      int tStep = bndry & TOPBNDRY ? 1 : nY;    // x -> x+h

      /* scorro sui punti della regione */
      for (int i = iInf; i <= iSup; ++i) {
        for (int j = jInf; j <= jSup; ++j) {
          int k = j - jInf;
          sum += std::abs(pde.residual(
              x[i], y[j], sol[i * nY + j], pde.rhs[i * nY + j],
              pde.neum & bndry, b[k], t[k], l[k], r[k], h));
        }

        /* muovo i puntatori da x a x+h */
        b += bStep;
        t += tStep;
        l += nY;
        r += nY;
      }
    }
  }

  return sum;
}

void PDESolver::init_point_bndry(int i, int j, const Real *sol,
                                 const PDE &pde,
                                 const Real *&b,
                                 const Real *&t,
                                 const Real *&l,
                                 const Real *&r) const {
  if (j == 0) { // siamo accanto al bordo inferiore
    b = pde.bottom + i;
  } else {
    b = sol + i * nY + j - 1;
  }

  if (j == nY - 1) { // siamo accanto al bordo superiore
    t = pde.top + i;
  } else {
    t = sol + i * nY + j + 1;
  }

  if (i == 0) { // siamo accanto al bordo sinistro
    l = pde.left + j;
  } else {
    l = sol + (i - 1) * nY + j;
  }

  if (i == nX - 1) { // siamo accanto al bordo destro
    r = pde.right + j;
  } else {
    r = sol + (i + 1) * nY + j;
  }
}
