#ifndef PDE_HPP
#define PDE_HPP

#include "BndryProp.hpp"
#include "Real.hpp"

struct PDE
{
  /**
   * Array-like container to store the value of the right-hand side of the PDE
   * at each point in lattice.
   */
  Real * rhs;

  /**
   * @brief Solves a 1-point PDE problem.
   * Returns the value solving the PDE at point (x,y) for some conditions on
   * the adjacent points and the given value of the right hand side.
   * @param neum Bitfield to flag the sides of the boundary where a Neumann
   * condition is set.
   * @param x The x coordinate.
   * @param y The y coordinate.
   * @param bottom The condition at (x, y-h).
   * @param top The condition at (x, y+h).
   * @param left The condition at (x-h, y).
   * @param right The condition at (x+h, y).
   * @param rhs The PDE right-hand side as evaluated at (x,y).
   */
  Real (*next_value)(Real x,
                     Real y,
                     Real rhs,
                     BndryProp neum,
                     Real bottom,
                     Real top,
                     Real left,
                     Real right,
                     Real h);

  /**
   * @brief Returns the residual at (x,y).
   * @param neum Bitfield to flag the sides of the boundary where a Neumann
   * condition is set.
   * @param x The x coordinate.
   * @param y The y coordinate.
   * @param middle The value at (x,y).
   * @param bottom The value at (x, y-h).
   * @param top The value at (x, y+h).
   * @param left The value at (x-h, y).
   * @param right The value at (x+h, y).
   * @param rhs Right-hand side of the PDE evaluated at (x,y).
   */
  Real (*residual)(Real x,
                   Real y,
                   Real middle,
                   Real rhs,
                   BndryProp neum,
                   Real bottom,
                   Real top,
                   Real left,
                   Real right,
                   Real h);

  /**
   * Specifies the sides of the boundary for which Neumann boundary conditions
   * are provided.
   */
  BndryProp neum;

  /**
   * Pointer to an array of nX elements containing the conditions for the top
   * side (where y is maximum) of the boundary, ordered by the value of x
   * they're associated with, from minimum to maximum.
   */
  Real * top;

  /**
   * Pointer to an array of nY elements containing the conditions for the left
   * side (where x is minimum) of the boundary, ordered by the value of y
   * they're associated with, from minimum to maximum.
   */
  Real * left;

  /**
   * Pointer to an array of nY elements containing the conditions for the right
   * side (where x is maximum) of the boundary, ordered by the value of y
   * they're associated with, from minimum to maximum.
   */
  Real * right;

  /**
   * Pointer to an array of nX elements containing the conditions for the
   * bottom side (where y is minimum) of the boundary, ordered by the value of
   * x they're associated with, from minimum to maximum.
   */
  Real * bottom;
};

#endif
