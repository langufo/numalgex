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

Real pot_cyl(Real r2, Real a2) {
  if (r2 <= a2) {
    return 0.25 * (a2 - r2);
  } else {
    return 0.25 * a2 * std::log(a2 / r2);
  }
}

/*
 * Argomenti da passare all'eseguibile:
 *  > dimensione della matrice
 *  > numero di iterazioni
 *  > algoritmo risolutivo: 0->SOR, 1->Gauss-Seidel, 2->Jacobi
 *  > correggere discontinuità? 0->no, altro->sì
 *  > prefisso per il nome dei file di output
 */
int main(int argc, char *argv[]) {
  using namespace std;

  long n = strtol(argv[1], nullptr,
                  0); // dimensione matrice
  long iter = strtol(argv[2], nullptr, 0); // iterazioni
  long algo = strtol(argv[3], nullptr, 0); // scelta algoritmo
  long correc = strtol(argv[4], nullptr, 0); // correzioni

  string prefix(argv[5]); // prefisso per il nome dei file
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
  vector<Real> an(n);      // soluzione esatta del problema

  Real h = static_cast<Real>(1) / (n + 1); // passo reticolo
  long a = n / 2 + 1; // raggio del cilindro in unità di h

  /* preparo condizione al bordo e soluzione esatta */
  for (long i = 0; i < n; ++i) {
    Real r2 = (i + 1) * (i + 1) * h * h;
    an[i] = pot_cyl(r2, a * a * h * h);
    nonzero[i] = pot_cyl(1, a * a * h * h);
  }

  /* popolo la matrice di rhs */
  Matrix<Real> m(n, rhs.data());
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i + 1 < a) {
        m[i][j] = -1;
      } else if (i + 1 > a) {
        m[i][j] = 0;
      } else { // i + 1 == a
        if (correc) {
          m[i][j] = -0.5; // rhs corretto per discontinuità
        } else {
          m[i][j] = -1;
        }
      }
    }
  }

  /* assemblo la definizione del problema */
  PDE pde;
  pde.rhs = rhs.data();
  pde.next_value = CylinderPDE::next_value;
  pde.residual = CylinderPDE::residual;
  pde.neum = BOTTOMBNDRY | LEFTBNDRY | TOPBNDRY;
  pde.bottom = zero.data();
  pde.top = zero.data();
  pde.left = zero.data();
  pde.right = nonzero.data();

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
      for (long i = 0; i < n; ++i) {
        Real x = (i + 1) * h;
        for (long j = 0; j < n; ++j) {
          Real y = (j + 1) * h;
          Real e = abs(m[i][j] - an[i]); // errore
          solFile << x << "\t" << y << "\t" << m[i][j] << "\n";
          errFile << x << "\t" << y << "\t" << e << "\n";
        }
        solFile << "\n";
        errFile << "\n";
      }
      solFile << endl;
      errFile << endl;
    }
  }

  delete solver;
}
