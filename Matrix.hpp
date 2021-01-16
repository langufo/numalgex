#ifndef MATRIX_HPP
#define MATRIX_HPP

/**
 * Interfaccia per accedere a un array con la sintassi di una
 * matrice
 */
template <typename T> class Matrix {
public:
  /**
   * @param cols Numero di colonne della matrice
   * @param data Puntatore al primo elemento
   */
  Matrix(int cols, T *data = nullptr) : n(cols), p(data) {}

  /**
   * @param i Riga cui si accede
   * @return Un puntatore al primo elemento della riga
   */
  T *operator[](int i) const { return p + i * n; }

  /**
   * Imposta l'indirizzo di memoria del primo elemento della
   * matrice.
   * @param data Puntatore al primo elemento
   */
  void set_first_elem(T *data) { p = data; }

private:
  /**
   * Numero di colonne della matrice
   */
  int n;

  /**
   * Primo elemento della matrice
   */
  T *p;
};

#endif
