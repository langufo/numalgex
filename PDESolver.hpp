#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include <vector>

#include "Matrix.hpp"

class PDESolver
{
public:
  /**
   * @brief Refine an approximated solution.
   * @param sol Pointer to the start of a contiguous region of memory holding
   * the initial guess for the solution to be overwritten. The values are
   * ordered by x first, and then by y.
   * @param resAsErr If true, the error reported by this method is the sum of
   * (the absolute values of) the residuals; otherwise, it's the sum of (the
   * absolute values of) the corrections applied at each point in the lattice.
   * @return The error of the refined solution as per resAsErr.
   */
  virtual double iter(double *sol, bool resAsErr) = 0;

  typedef unsigned BndryLayout;
  static const BndryLayout NOBNDRY = 0;
  static const BndryLayout TOPBNDRY = 1;
  static const BndryLayout LEFTBNDRY = 2;
  static const BndryLayout RIGHTBNDRY = 4;
  static const BndryLayout BOTTOMBNDRY = 8;

protected:
  /**
   * @param rhs Pointer to a function of the coordinates which acts as the right
   * hand side of the PDE.
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
   * @param left Pointer to an array of nY elements containing the conditions on
   * the left side of the boundary.
   * @param right Like @p left, but for the right side.
   */
  PDESolver(double (*rhs)(double, double), BndryLayout neumLayout, double x0,
            double y0, double h, int nX, int nY, const double *bottom,
            const double *top, const double *left, const double *right);

  virtual double next_value(BndryLayout neum, double x, double y, double bottom,
                            double top, double left, double right,
                            double rhs) const = 0;

  virtual double residual(BndryLayout neum, double x, double y, double middle,
                          double bottom, double top, double left, double right,
                          double rhs) const = 0;

  double line_error(int i, int jFirst, int jLast, BndryLayout neum,
                    const double *bottom, const double *top, const double *left,
                    const double *right) const;

  double abs_res_sum() const;

  /**
   * Distance between adjacent points in the lattice.
   */
  double h;

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

  std::vector<double> x, y;

  /**
   * Specifies the sides of the boundary for which Neumann boundary conditions
   * are provided.
   */
  BndryLayout neumLayout;

  /**
   * Pointer to an array of nX elements containing the conditions for the bottom
   * side (where y is minimum) of the boundary, ordered by the value of x
   * they're associated with, from minimum to maximum.
   */
  const double *bottom;

  /**
   * Pointer to an array of nX elements containing the conditions for the top
   * side (where y is maximum) of the boundary, ordered by the value of x
   * they're associated with, from minimum to maximum.
   */
  const double *top;

  /**
   * Pointer to an array of nY elements containing the conditions for the left
   * side (where x is minimum) of the boundary, ordered by the value of y
   * they're associated with, from minimum to maximum.
   */
  const double *left;

  /**
   * Pointer to an array of nY elements containing the conditions for the right
   * side (where x is maximum) of the boundary, ordered by the value of y
   * they're associated with, from minimum to maximum.
   */
  const double *right;

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
  BndryLayout bndryX[3];

  /**
   * Specifies for each of the regions along the y axis whether they border the
   * top or bottom boundary.
   */
  BndryLayout bndryY[3];

  /**
   * The first and last lattice row index for each of the regions along the x
   * axis.
   */
  int limX[3][2];

  /**
   * The first and last lattice column index for each of the regions along the y
   * axis.
   */
  int limY[3][2];

  Matrix<double> ms;

private:
  void cache_rhs(double (*rhs)(double, double));

  std::vector<double> r;

protected:
  Matrix<double> mr;
};

#endif
