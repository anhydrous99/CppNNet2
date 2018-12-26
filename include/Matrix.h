//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_MATRIX_H
#define CPPNNET2_MATRIX_H

#include <cstddef>
#include <memory>

namespace CppNNet2 {
  template<class T>
  class Matrix {
  public:
    // types:
    typedef T value_type;

  private:
    size_t _rows, _cols, _size;
    // Data pointers
    value_type *begin;
    value_type *end;

  public:
    // construct/destroy
    Matrix() : begin(0), end(0) { _rows = _cols = _size = 0; }

    inline explicit Matrix(size_t n);
    Matrix(const value_type &x, size_t n);
    Matrix(const value_type *px, size_t n);
    Matrix(const Matrix &m);
    Matrix(Matrix &&m) noexcept;
    Matrix(size_t n, size_t m);
    Matrix(const value_type &x, size_t n, size_t m);
    Matrix(const value_type &px, size_t n, size_t m);
    ~Matrix();

    // assignment
    Matrix &operator=(const Matrix &m);
    Matrix &operator=(Matrix &&m) noexcept;
    Matrix &operator=(const value_type &x);

    // element access
    const value_type &operator[](size_t i) const;
    value_type &operator[](size_t i);
    const value_type &operator()(size_t i, size_t j) const;
    value_type &operator()(size_t i, size_t j);
    const value_type cat(size_t i) const;
    value_type &at(size_t i);
    const value_type &cat(size_t i, size_t j) const;
    value_type &at(size_t i, size_t j);

    // unary operators:
    Matrix operator+() const;
    Matrix operator-() const;
    Matrix operator~() const;
    Matrix<bool> operator!() const;

    // computed assignment:
    Matrix &operator*=(const value_type &x);
    Matrix &operator/=(const value_type &x);
    Matrix &operator%=(const value_type &x);
    Matrix &operator+=(const value_type &x);
    Matrix &operator-=(const value_type &x);
    Matrix &operator^=(const value_type &x);
    Matrix &operator&=(const value_type &x);
    Matrix &operator|=(const value_type &x);
    Matrix &operator<<=(const value_type &x);
    Matrix &operator>>=(const value_type &x);

    Matrix &operator*=(const Matrix &m);
    Matrix &operator/=(const Matrix &m);
    Matrix &operator%=(const Matrix &m);
    Matrix &operator+=(const Matrix &m);
    Matrix &operator-=(const Matrix &m);
    Matrix &operator^=(const Matrix &m);
    Matrix &operator|=(const Matrix &m);
    Matrix &operator&=(const Matrix &m);
    Matrix &operator<<=(const Matrix &m);
    Matrix &operator>>=(const Matrix &m);

    void swap(Matrix &v) noexcept;

    constexpr size_t size();
    constexpr size_t rows();
    constexpr size_t cols();

    void clear(size_t capacity);
    void resize(size_t n, value_type x = value_type());
    void resize(size_t n, size_t m, value_type x = value_type());
  };

  template<class T>
  void swap(Matrix<T> &x, Matrix<T> y);

  template<class T>
  Matrix<T> operator*(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator*(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator*(const T &y, const Matrix<T> &x);

  template<class T>
  Matrix<T> operator/(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator/(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator/(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> operator%(const Matrix<T> &x, const Matrix<T> y);
  template<class T>
  Matrix<T> operator%(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator%(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> operator+(const Matrix<T> &x, const Matrix<T> y);
  template<class T>
  Matrix<T> operator+(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator+(const T &y, const Matrix<T> &x);

  template<class T>
  Matrix<T> operator-(const Matrix<T> &x, const Matrix<T> y);
  template<class T>
  Matrix<T> operator-(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator-(const T &y, const Matrix<T> &x);

  template<class T>
  Matrix<T> operator^(const Matrix<T> &x, const Matrix<T> y);
  template<class T>
  Matrix<T> operator^(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator^(const T &y, const Matrix<T> &y);

  template<class T>
  Matrix<T> operator&(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator&(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator&(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> operator|(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator|(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator|(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> operator<<(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator<<(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator<<(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> operator>>(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator>>(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator>>(const T &x, Matrix<T> &y);

  template<class T>
  Matrix<bool> operator&&(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator&&(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator&&(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<bool> operator||(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator||(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator||(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<bool> operator==(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator==(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator==(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<bool> operator!=(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator!=(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator!=(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<bool> operator<(const Matrix<T> &x, const Matrix<T> y);
  template<class T>
  Matrix<bool> operator<(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator<(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<bool> operator>(const Matrix<T> x, const Matrix<T> y);
  template<class T>
  Matrix<bool> operator>(const Matrix<T> x, const T &y);
  template<class T>
  Matrix<bool> operator>(const T &x, const Matrix<T> y);

  template<class T>
  Matrix<bool> operator<=(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator<=(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator<=(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<bool> operator>=(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator>=(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator>=(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> abs(const Matrix<T> &x);
  template<class T>
  Matrix<T> acos(const Matrix<T> &x);
  template<class T>
  Matrix<T> asin(const Matrix<T> &x);

  template<class T>
  Matrix<T> atan2(const Matrix<T> &x, const Matrix<T> y);
  template<class T>
  Matrix<T> atan2(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> atan2(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> cos(const Matrix<T> &x);
  template<class T>
  Matrix<T> cosh(const Matrix<T> &x);
  template<class T>
  Matrix<T> exp(const Matrix<T> &x);
  template<class T>
  Matrix<T> log(const Matrix<T> &x);
  template<class T>
  Matrix<T> log10(const Matrix<T> &x);

  template<class T>
  Matrix<T> pow(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> pow(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> pow(const T &x, const Matrix<T> &y);

  template<class T>
  Matrix<T> sin(const Matrix<T> &x);
  template<class T>
  Matrix<T> sinh(const Matrix<T> &x);
  template<class T>
  Matrix<T> sqrt(const Matrix<T> &x);
  template<class T>
  Matrix<T> tan(const Matrix<T> &x);
  template<class T>
  Matrix<T> tanh(const Matrix<T> &x);

// end

  template<class T>
  inline Matrix<T>::Matrix(size_t n) : begin(0), end(0) {
    _cols = _size = n;
    _rows = 1;
    if (n) {
      begin = end = static_cast<value_type *>(malloc(sizeof(value_type) * n));
      try {
        for (size_t n_left = n; n_left; --n_left, ++end)
          new(end) value_type();
      } catch (...) {
        clear(n);
        throw;
      }
    }
  }

  template<class T>
  inline Matrix<T>::Matrix(const value_type &x, size_t n) : begin(0), end(0) {
    resize(n, x);
  }

  template<class T>
  Matrix<T>::Matrix(const value_type *p, size_t n) {
  }
}

#endif //CPPNNET2_MATRIX_H
