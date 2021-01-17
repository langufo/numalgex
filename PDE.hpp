#ifndef PDE_HPP
#define PDE_HPP

#include "BndryProp.hpp"
#include "Real.hpp"

struct PDE {
  /**
   * Puntatore a un array che contiene, per ciascun punto del
   * reticolo, ordinati prima per x e poi per y, i right-hand
   * side dell'equazione
   */
  Real *rhs;

  /**
   * Determina il valore della funzione in un punto.
   * @param x Prima coordinata
   * @param y Seconda coordinata
   * @param rhs Right-hand side dell'equazione, valutato nel
   * punto (x,y)
   * @param neum Indica i lati per i quali è fornita la
   * derivata invece che valore della funzione
   * @param bottom Valore per il punto (x,y-h)
   * @param top Valore per il punto (x,y+h)
   * @param left Valore per il punto (x-h,y)
   * @param right Valore per il punto (x+h,y)
   * @param h Passo del reticolo
   * @return Il valore della funzione determinato nel punto
   * (x,y).
   */
  Real (*next_value)(Real x, Real y, Real rhs, BndryProp neum,
                     Real bottom, Real top, Real left,
                     Real right, Real h);

  /**
   * Calcola il residuo in un punto.
   * @param x Prima coordinata
   * @param y Seconda coordinata
   * @param middle Valore della funzione nel punto (x,y)
   * @param rhs Right-hand side dell'equazione, valutato nel
   * punto (x,y)
   * @param neum Indica i lati per i quali è fornita la
   * derivata al posto del valore della funzione
   * @param bottom Valore per il punto (x,y-h)
   * @param top Valore per il punto (x,y+h)
   * @param left Valore per il punto (x-h,y)
   * @param right Valore per il punto (x+h,y)
   * @param h Passo del reticolo
   * @return Residuo per il punto (x,y).
   */
  Real (*residual)(Real x, Real y, Real middle, Real rhs,
                   BndryProp neum, Real bottom, Real top,
                   Real left, Real right, Real h);

  /**
   * Identifica i lati per i quali sono fornite condizioni
   * al contorno di Neumann.
   */
  BndryProp neum;

  /**
   * Puntatore a un array che ospita la parte di condizioni al
   * contorno relativa al bordo inferiore (y minimo)
   */
  Real *bottom;

  /**
   * Puntatore a un array che ospita la parte di condizioni al
   * contorno relativa al bordo superiore (y massimo)
   */
  Real *top;

  /**
   * Puntatore a un array che ospita la parte di condizioni al
   * contorno relativa al bordo sinistro (x minimo)
   */
  Real *left;

  /**
   * Puntatore a un array che ospita la parte di condizioni al
   * contorno relativa al bordo destro (x massimo)
   */
  Real *right;
};

#endif
