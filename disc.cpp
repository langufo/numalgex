#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "BndryProp.hpp"
#include "CylinderPDE.hpp"
#include "GSeidelSolver.hpp"
#include "JacobiSolver.hpp"
#include "Matrix.hpp"
#include "PDE.hpp"
#include "PDESolver.hpp"
#include "Real.hpp"
#include "SORSolver.hpp"

Real pot_axis(Real z, Real r2, Real height) {
  z = std::abs(z);
  Real a = z - height / 2;
  Real b = z + height / 2;
  Real p = std::sqrt(r2 + a * a);
  Real q = std::sqrt(r2 + b * b);
  if (a < 0) {
    return 0.25 *
           (b * q - a * p - 2 * z * z - height * height / 2 -
            r2 * std::log(r2 / ((p - a) * (q + b))));
  } else {
    return 0.25 * (height * (2 * z * z / (p + q) +
                             (p + q) / 2 - 2 * z) -
                   r2 * std::log((p + a) / (q + b)));
  }
}

/*
 * Argomenti da passare all'eseguibile:
 *  > dimensione della matrice
 *  > numero di iterazioni
 *  > algoritmo risolutivo: 0->SOR, 1->Gauss-Seidel, 2->Jacobi
 *  > prefisso per il nome dei file di output
 */
int main(int argc, char *argv[]) {
  using namespace std;

  long n = strtol(argv[1], nullptr, 0); // dimensione matrice
  long iter = strtol(argv[2], nullptr, 0); // iterazioni
  long algo = strtol(argv[3], nullptr, 0); // scelta algoritmo

  string prefix(argv[4]); // prefisso per il nome dei file
  ofstream iterFile(prefix + "iter.txt"); // info iterazioni
  ofstream solFile(prefix + "sol.txt");   // soluzioni complete
  ofstream axisFile(prefix + "axis.txt"); // soluzioni su asse
  ofstream errFile(prefix + "err.txt");   // errori sull'asse
  scientific(iterFile);
  scientific(solFile);
  scientific(axisFile);
  scientific(errFile);

  vector<Real> zero(n);    // condizione al bordo nulla
  vector<Real> nonzero(n); // condizione al bordo non nulla
  vector<Real> sol(n * n); // memoria per la soluzione
  vector<Real> rhs(n * n); // memoria per rhs
  vector<Real> an(n);      // soluzione esatta sull'asse

  Real h = static_cast<Real>(1) / (n + 1); // passo reticolo
  long a = n / 8 + 1;           // raggio in unità di h
  long halfHeight = n / 16 + 1; // mezza altezza in unità di h
  Real h2 = h * h;
  long a2 = a * a;

  /* preparo condizione al bordo e soluzione esatta */
  for (long j = 0; j < n; ++j) {
    an[j] = pot_axis((j + 1) * h, a2 * h2, 2 * halfHeight * h);
    nonzero[j] = 0.5 * a2 * halfHeight * h2 /
                 sqrt((j + 1) * (j + 1) + (n + 1) * (n + 1));
  }

  /* popolo la matrice di rhs */
  Matrix<Real> m(n, rhs.data());
  for (long i = 0; i < n; ++i) {
    for (long j = 0; j < n; ++j) {
      if (i + 1 > a || j + 1 > halfHeight) {
        m[i][j] = 0;
      } else {
        m[i][j] = -1;
      }
      if (i + 1 == a) {
        m[i][j] /= 2; // correzione rhs per faccia esterna
      }
      if (j + 1 == halfHeight) {
        m[i][j] /= 2; // correzione rhs per faccia superiore
      }
    }
  }

  /* assemblo la definizione del problema */
  PDE pde;
  pde.rhs = rhs.data();
  pde.next_value = CylinderPDE::next_value;
  pde.residual = CylinderPDE::residual;
  pde.neum = BOTTOMBNDRY | LEFTBNDRY;
  pde.bottom = zero.data();
  pde.left = zero.data();
  pde.right = nonzero.data();
  pde.top = nonzero.data();

  PDESolver *solver = nullptr;
  switch (algo) {
  case 0: {
    Real w = 2 / (1 + sqrt(4.45) / (n + 2));
    solver = new SORSolver(h, h, h, n, n, w);
    break;
  }
  case 1:
    solver = new GSeidelSolver(h, h, h, n, n);
    break;
  case 2:
    solver = new JacobiSolver(h, h, h, n, n);
    break;
  }

  m.set_first_elem(sol.data()); // cambio di matrice!
  Real thrErr = 1; // soglia per stampa (errore iterazione)
  Real thrRes = 1; // altra soglia per stampa (residuo)
  for (int k = 0; k < iter; ++k) {
    Real err = solver->iter(sol.data(), pde);
    Real res = solver->abs_res_sum(sol.data(), pde);
    iterFile << k + 1 << "\t" << err << "\t" << res << "\n";

    if (err != 0 && (err < thrErr || res < thrRes)) { // stampa
      /* aggiorno la soglia raggiunta decimandola */
      if (err < thrErr) {
        thrErr /= 10;
      }
      if (res < thrRes) {
        thrRes /= 10;
      }

      /* intestazione */
      solFile << "# " << k + 1 << " " << err << " " << res
              << "\n";
      axisFile << "# " << k + 1 << " " << err << " " << res
               << "\n";
      errFile << "# " << k + 1 << " " << err << " " << res
              << "\n";
      /* soluzione ed errore */
      for (long j = 0; j < n; ++j) {
        Real y = (j + 1) * h;
        solFile << "\n";
        for (long i = 0; i < n; ++i) {
          Real x = (i + 1) * h;
          solFile << x << "\t" << y << "\t" << m[i][j] << "\n";
        }
        axisFile << y << "\t" << m[0][j] << "\n";
        Real e = abs(m[0][j] - an[j]); // errore
        errFile << y << "\t" << e << "\n";
      }
      solFile << "\n\n";
      axisFile << "\n\n";
      errFile << "\n\n";
      flush(solFile);
      flush(axisFile);
      flush(errFile);
    }
  }

  delete solver;
}
