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

Real pot_sphere(Real r2, Real a2) {
  if (r2 <= a2) {
    return 0.5 * (a2 - r2 / 3);
  } else {
    return a2 * std::sqrt(a2 / r2) / 3;
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
  ofstream solFile(prefix + "sol.txt");   // soluzioni trovate
  ofstream errFile(prefix + "err.txt");   // errori soluzioni
  scientific(iterFile);
  scientific(solFile);
  scientific(errFile);

  vector<Real> zero(n);    // condizione al bordo nulla
  vector<Real> nonzero(n); // condizione al bordo non nulla
  vector<Real> sol(n * n); // memoria per la soluzione
  vector<Real> rhs(n * n); // memoria per rhs

  Real h = static_cast<Real>(1) / (n + 1); // passo reticolo
  long a = n / 2 + 1; // raggio della sfera in unità di h
  Real h2 = h * h;
  long a2 = a * a;

  /* preparo la condizione al bordo non nulla */
  for (long i = 0; i < n; ++i) {
    nonzero[i] = pot_sphere(
        ((i + 1) * (i + 1) + (n + 1) * (n + 1)) * h2, a2 * h2);
  }

  /* popolo la matrice di rhs */
  Matrix<Real> m(n, rhs.data());
  for (long i = 0; i < n; ++i) {
    for (long j = 0; j < n; ++j) {
      if ((i + 1) * (i + 1) + (j + 1) * (j + 1) <= a * a) {
        m[i][j] = -1;
      } else {
        m[i][j] = 0;
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
      errFile << "# " << k + 1 << " " << err << " " << res
              << "\n";
      /* soluzione ed errore */
      for (int j = 0; j < n; ++j) {
        solFile << "\n";
        errFile << "\n";
        Real y = (j + 1) * h;
        for (int i = 0; i < n; ++i) {
          Real x = (i + 1) * h;
          Real e = abs(m[i][j] -
                       pot_sphere(x * x + y * y, a2 * h2));
          solFile << x << "\t" << y << "\t" << m[i][j] << "\n";
          errFile << x << "\t" << y << "\t" << e << "\n";
        }
      }
      solFile << "\n\n";
      errFile << "\n\n";
      flush(solFile);
      flush(errFile);
    }
  }

  delete solver;
}
