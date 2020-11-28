#ifndef MATRIX_HPP
#define MATRIX_HPP

template<typename T>
class Matrix
{
public:
  Matrix(int cols, T * data = nullptr)
    : n(cols)
    , p(data)
  {}

  T * operator[](int i) const { return p + i * n; }

  void set_first_elem(T * data) { p = data; }

private:
  int n;
  T * p;
};

#endif
