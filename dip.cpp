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

Real pot_axis(Real z, Real a) {
  if (z == 0) {
    return 0;
  }
  bool neg = z < 0;
  Real z2 = z * z;
  Real a2 = a * a; // quadrato del raggio della sfera
  Real r2 = a2 + z2;
  Real r = std::sqrt(r2); // norma della distanza
  z = std::abs(z);
  Real v;
  if (z > a) {
    v = (r2 * r / z - z2) / 3 - a2 / 2;
  } else {
    v = (r2 * r - a2 * a) / (3 * z) - z2 / 2;
  }
  return neg ? -v : v;
}

/*
 * Lista degli argomenti:
 *  - dimensione della matrice
 *  - numero di iterazioni
 *  - algoritmo risolutivo: 0 -> SOR, 1 -> Gauss-Seidel, 2 ->
 * Jacobi
 *  - prefisso per il nome dei file di output
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
  vector<Real> r(n);       // condizione al bordo destro
  vector<Real> t(n);       // condizione al bordo superiore
  vector<Real> sol(n * n); // memoria per la soluzione
  vector<Real> rhs(n * n); // memoria per rhs
  vector<Real> an(n);      // soluzione esatta sull'asse

  Real h = static_cast<Real>(1) / (n + 1); // passo reticolo
  long a = n / 8 + 1; // raggio della sfera in unità di h

  /* preparo condizioni al bordo e soluzione esatta */
  for (int i = 0; i < n; ++i) {
    long d2 = (i + 1) * (i + 1) + (n + 1) * (n + 1);
    Real coeff = a * a * a * a * h * h / 8;
    an[i] = pot_axis((i + 1) * h, a * h);
    r[i] = coeff * (i + 1) / (d2 * sqrt(d2));
    t[i] = coeff * (n + 1) / (d2 * sqrt(d2));
  }

  /* popolo la matrice di rhs */
  Matrix<Real> m(n, rhs.data());
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      long r = i + 1;
      long z = j + 1;
      if (r * r + z * z <= a * a) {
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
  pde.neum = LEFTBNDRY;
  pde.bottom = zero.data();
  pde.left = zero.data();
  pde.right = r.data();
  pde.top = t.data();

  PDESolver *solver = nullptr;
  switch (algo) {
  case 0: {
    Real w = 2 / (1 + std::sqrt(4.45) / (n + 2));
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
          solFile << (i + 1) * h << "\t" << y << "\t"
                  << m[i][j] << "\n";
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
