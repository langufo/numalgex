#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include <vector>

#include "BndryProp.hpp"
#include "PDE.hpp"
#include "Real.hpp"

class PDESolver
{
public:
  virtual ~PDESolver() {}

  /**
   * @brief Refine an approximated solution.
   * @param sol Pointer to an array holding the initial guess for the solution
   * to be overwritten. The values are ordered by x first, and then by y.
   * @param resAsErr If true, the error reported by this method is the sum of
   * (the absolute values of) the residuals; otherwise, it's the sum of (the
   * absolute values of) the corrections applied at each point in the lattice.
   * @return The sum of the absolute values of the corrections applied at each
   * point in the lattice.
   */
  virtual Real iter(Real * sol, const PDE & pde) = 0;

  /**
   * @return The sum over all the points in the lattice of the absolute values
   * of the residuals.
   */
  Real abs_res_sum(Real * sol, const PDE & pde) const;

protected:
  /**
   * @param rhs Pointer to an array containing the values of the right-hand
   * side of the PDE at the inner points of the lattice, ordered by x first,
   * and then by y.
   * @param neumLayout Specifies for which sides of the boundary Neumann
   * boundary conditions are provided.
   * @param x0 The minimum value of x for an inner point.
   * @param y0 The minimum value of y for an inner point.
   * @param h The minimum distance between two points along each axis.
   * @param nX The number of inner points along the x axis.
   * @param nY The number of inner points along the y axis.
   * @param bottom Pointer to an array of nX elements containing the conditions
   * on the bottom side of the boundary.
   * @param top Like @p bottom, but for the top side.
   * @param left Pointer to an array of nY elements containing the conditions
   * on the left side of the boundary.
   * @param right Like @p left, but for the right side.
   */
  PDESolver(Real x0, Real y0, Real h, int nX, int nY);

  /**
   * @param i Index of the desired x coordinate inside vector x.
   * @param j Index of the desired y coordinate inside vector y.
   * @param bndry Bitfield to flag the sides of the boundary the point at (x,y)
   * is a neighbour of.
   * @param b Reference to a pointer to point to the value at (x, y-h).
   * @param t Reference to a pointer to point to the value at (x, y+h).
   * @param l Reference to a pointer to point to the value at (x-h, y).
   * @param r Reference to a pointer to point to the value at (x+h, y).
   * @param bStep Offset that makes b point to the value at (x+h, y-h).
   * @param tStep Offset that makes t point to the value at (x+h, y+h).
   */
  void init_point_bndry(int i,
                        int j,
                        const Real * sol,
                        const PDE & pde,
                        const Real *& b,
                        const Real *& t,
                        const Real *& l,
                        const Real *& r) const;

  /**
   * Distance between adjacent points in the lattice.
   */
  Real h;

  /**
   * Number of inner points along the x axis; alternatively, number of rows in
   * the lattice used to represent the solution to the PDE.
   */
  int nX;

  /**
   * Number of inner points along the y axis; alternatively, number of columns
   * in the lattice used to represent the solution to the PDE.
   */
  int nY;

#ifndef ON_THE_FLY
  std::vector<Real> x, y;
#else
  class OTF
  {
  public:
    OTF(Real z0, Real h)
      : z0(z0)
      , h(h)
    {}
    Real operator[](int k) const { return z0 + k * h; }
    Real z0, h;
  } x, y;
#endif

  /**
   * Number of regions the lattice can be divided into along the x axis. They
   * shall not be more than 3.
   */
  int nRegX;

  /**
   * Number of regions the lattice can be divided into along the y axis. They
   * shall not be more than 3.
   */
  int nRegY;

  /**
   * Specifies for each of the regions along the x axis whether they border the
   * left or right boundary.
   */
  BndryProp bndryX[3];

  /**
   * Specifies for each of the regions along the y axis whether they border the
   * top or bottom boundary.
   */
  BndryProp bndryY[3];

  /**
   * The first and last lattice row index for each of the regions along the x
   * axis.
   */
  int limX[3][2];

  /**
   * The first and last lattice column index for each of the regions along the
   * y axis.
   */
  int limY[3][2];
};

#endif
