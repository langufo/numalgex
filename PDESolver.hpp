#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include <vector>

#include "Matrix.hpp"

class PDESolver
{
public:
  /**
   * @param sol Pointer to the start of a contiguous region of memory holding
   * an initial solution to the problem used to derive another that then
   * replaces the previous one.
   * @param resAsErr If true, the error reported by this method is the sum of
   * (the absolute values of) the residuals; otherwise, it's the sum of (the
   * absolute values of) the corrections applied at each point in the lattice.
   */
  virtual double iter(double* sol, bool resAsErr) = 0;

  typedef unsigned BndryLayout;
  static const BndryLayout NOBNDRY = 0;
  static const BndryLayout TOPBNDRY = 1;
  static const BndryLayout LEFTBNDRY = 2;
  static const BndryLayout RIGHTBNDRY = 4;
  static const BndryLayout BOTTOMBNDRY = 8;

protected:
  PDESolver(double (*rhs)(double, double),
            BndryLayout neumLayout,
            double x0,
            double y0,
            double h,
            int nX,
            int nY,
            const double* bottom,
            const double* top,
            const double* left,
            const double* right);

  virtual double next_value(BndryLayout neum,
                            double x,
                            double y,
                            double bottom,
                            double top,
                            double left,
                            double right,
                            double rhs) = 0;
  virtual double residual(BndryLayout neum,
                          double x,
                          double y,
                          double middle,
                          double bottom,
                          double top,
                          double left,
                          double right,
                          double rhs) = 0;

  double line_error(int i,
                    int jFirst,
                    int jLast,
                    BndryLayout neum,
                    const double* bottom,
                    const double* top,
                    const double* left,
                    const double* right);

  double abs_res_sum();

  double (*rhs)(double, double);

  double h;
  int nX, nY;
  BndryLayout neumLayout;
  std::vector<double> x, y;

  const double* bottom;
  const double* top;
  const double* left;
  const double* right;

  int nRegX, nRegY;
  BndryLayout bndryX[3], bndryY[3];
  int limX[3][2], limY[3][2];

private:
  void cache_rhs();

  std::vector<double> r;

protected:
  Matrix<double> mr, ms;
};

#endif
