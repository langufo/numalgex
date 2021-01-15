#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include <vector>

#include "BndryProp.hpp"
#include "PDE.hpp"
#include "Real.hpp"

class PDESolver {
public:
  virtual ~PDESolver() {}

  /**
   * Raffina la soluzione di un problema.
   * @param sol Puntatore a un array contenente la soluzione di
   * partenza, che sarà sovrascritta
   * @param pde Definizione del problema
   * @return La somma dei valori assoluti delle correzioni
   * apportate in ciascun punto ("errore dell'iterazione").
   */
  virtual Real iter(Real *sol, const PDE &pde) = 0;

  /**
   * @param sol Puntatore a un array contenente la soluzione
   * del problema
   * @return La somma estesa a tutti i punti del reticolo dei
   * valori assoluti dei residui.
   */
  Real abs_res_sum(Real *sol, const PDE &pde) const;

protected:
  /**
   * @param x0 Valore minimo della prima coordinata di un punto
   * del reticolo
   * @param y0 Valore minimo della seconda coordinata per un
   * punto del reticolo
   * @param h Passo del reticolo
   * @param nX Numero di punti nella prima direzione
   * @param nY Numero di punti nella seconda direzione
   */
  PDESolver(Real x0, Real y0, Real h, int nX, int nY);

  /**
   * Metodo d'utilità che restituisce i puntatori agli elementi
   * adiacenti a una posizione del reticolo.
   * @param i Primo indice del punto centrale
   * @param j Secondo indice del punto centrale
   * @param sol Puntatore all'array contenente la soluzione
   * @param pde Definizione del problema
   * @param b Punta alla posizione sottostante
   * @param t Punta alla posizione soprastante
   * @param l Punta alla posizione a sinistra
   * @param r Punta alla posizione a destra
   */
  void init_point_bndry(int i, int j, const Real *sol,
                        const PDE &pde, const Real *&b,
                        const Real *&t, const Real *&l,
                        const Real *&r) const;

  /**
   * Passo del reticolo
   */
  Real h;

  /**
   * Numero di punti nella prima direzione
   */
  int nX;

  /**
   * Numero di punti nella seconda direzione
   */
  int nY;

  /**
   * Prime coordinate dei punti del reticolo
   */
  std::vector<Real> x;

  /**
   * Seconde coordinate dei punti del reticolo
   */
  std::vector<Real> y;

  /**
   * Numero di divisioni del reticolo nella prima direzione
   */
  int nRegX;

  /**
   * Numero di divisioni del reticolo nella seconda direzione
   */
  int nRegY;

  /**
   * Per ciascuna divisione nella prima direzione specifica se
   * essa è adiacente al bordo sinistro o destro
   */
  BndryProp bndryX[3];

  /**
   * Per ciascuna divisione nella seconda direzione specifica
   * se essa è adiacente al bordo inferiore o superiore
   */
  BndryProp bndryY[3];

  /**
   * Gli indici delle estremità di ciascuna suddivisione nella
   * prima direzione
   */
  int limX[3][2];

  /**
   * Gli indici delle estremità di ciascuna suddivisione nella
   * seconda direzione
   */
  int limY[3][2];
};

#endif
