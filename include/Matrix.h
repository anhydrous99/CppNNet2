//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_MATRIX_H
#define CPPNNET2_MATRIX_H

#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <memory>

#ifdef USES_OPENBLAS

#include <cblas.h>

#endif

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

    explicit Matrix(size_t n);
    Matrix(const value_type *px, size_t n);
    Matrix(const Matrix &m);
    Matrix(Matrix &&m) noexcept;
    Matrix(size_t n, size_t m);
    Matrix(const value_type &x, size_t n, size_t m);
    Matrix(const value_type *px, size_t n, size_t m);
    ~Matrix();

    // assignment
    Matrix &operator=(const Matrix &m);
    Matrix &operator=(Matrix &&m) noexcept;
    Matrix &operator=(const value_type &x);

    // element access
    const value_type &operator[](size_t i) const { return begin[i]; }

    value_type &operator[](size_t i) { return begin[i]; }

    const value_type &operator()(size_t i, size_t j) const { return begin[i * _cols + j]; }

    value_type &operator()(size_t i, size_t j) { return begin[i * _cols + j]; }

    const value_type cat(size_t i) const { return begin[i]; }

    value_type &at(size_t i) { return begin[i]; }

    const value_type &cat(size_t i, size_t j) const { return begin[i * _cols + j]; }

    value_type &at(size_t i, size_t j) { return begin[i * _cols + j]; }

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

    constexpr size_t size() const { return _size; }

    constexpr size_t rows() const { return _rows; }

    constexpr size_t cols() const { return _cols; }

    void clear(size_t capacity);
    void resize(size_t n, value_type x = value_type());
    void resize(size_t n, size_t m, value_type x = value_type());

    T *data() { return begin; }

    const T *data() const { return begin; }

  private:
    Matrix &assign_range(const value_type *f, const value_type *l);
  };

  template<class T>
  Matrix<T> operator*(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator*(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator*(const T &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator/(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator/(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator/(const T &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator%(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator%(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator%(const T &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator+(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator+(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator+(const T &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator-(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator-(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator-(const T &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator^(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> operator^(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<T> operator^(const T &x, const Matrix<T> &y);
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
  Matrix<bool> operator<(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator<(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator<(const T &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator>(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<bool> operator>(const Matrix<T> &x, const T &y);
  template<class T>
  Matrix<bool> operator>(const T &x, const Matrix<T> &y);
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
  Matrix<T> atan2(const Matrix<T> &x, const Matrix<T> &y);
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
  template<class T>
  bool all_equal(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> transpose(const Matrix<T> &x);
  template<class T>
  Matrix<T> kronecker_product(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> box_product(const Matrix<T> &x, const Matrix<T> &y);
  template<class T>
  Matrix<T> hadamard_product(const Matrix<T> &x, const Matrix<T> &y);

  enum CONCAT_TYPE {
    HORIZONTAL, VERTICAL
  };

  template<class T>
  Matrix<T> concat(CONCAT_TYPE contype, const Matrix<T> &x, const Matrix<T> &y);

// end

  template<class T>
  Matrix<T>::Matrix(size_t n) : begin(0), end(0) {
    _cols = _size = n;
    _rows = 1;
    if (n) {
      begin = end = static_cast<value_type *>(std::malloc(sizeof(value_type) * n));
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
  Matrix<T>::Matrix(const value_type *p, size_t n) {
    _cols = _size = n;
    _rows = 1;
    begin = p;
    end = p + n;
  }

  template<class T>
  Matrix<T>::Matrix(const Matrix<T> &m) : begin(nullptr), end(nullptr) {
    _rows = m.rows();
    _cols = m.cols();
    _size = m.size();
    if (m.size()) {
      begin = end = static_cast<value_type *>(std::malloc(sizeof(value_type) * m.size()));
      try {
        for (value_type *p = m.begin; p != m.end; ++end, ++p)
          new(end) value_type(*p);
      } catch (...) {
        clear(m.size());
        throw;
      }
    }
  }

  template<class T>
  Matrix<T>::Matrix(Matrix<T> &&m) noexcept : begin(m.begin), end(m.end) {
    _rows = m.rows();
    _cols = m.cols();
    _size = m.size();
    m.begin = m.end = nullptr;
  }

  template<class T>
  Matrix<T>::Matrix(size_t n, size_t m) : begin(nullptr), end(nullptr) {
    _rows = n;
    _cols = m;
    _size = n * m;
    if (n && m) {
      begin = end = static_cast<value_type *>(std::malloc(sizeof(value_type) * _size));
      try {
        for (size_t n_left = _size; n_left; --n_left, ++end)
          new(end) value_type();
      } catch (...) {
        clear(n);
        throw;
      }
    }
  }

  template<class T>
  Matrix<T>::Matrix(const value_type &x, size_t n, size_t m) : begin(0), end(0) {
    resize(n, m, x);
  }

  template<class T>
  Matrix<T>::Matrix(const value_type *px, size_t n, size_t m) {
    _rows = n;
    _cols = m;
    _size = n * m;
    begin = px;
    end = px + _size;
  }

  template<class T>
  Matrix<T>::~Matrix() {
    clear(_size);
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator=(const Matrix<T> &m) {
    if (this != &m)
      return assign_range(m.begin, m.end);
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator=(Matrix<T> &&m) noexcept {
    clear(size());
    begin = m.begin;
    end = m.end;
    m.begin = nullptr;
    m.end = nullptr;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator=(const value_type &x) {
    std::fill(begin, end, x);
    return *this;
  }

  template<class T>
  Matrix<T> Matrix<T>::operator+() const {
    Matrix<value_type> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(std::malloc(sizeof(value_type) * _size));
      for (const value_type *p = begin; n; ++r.end, ++p, --n)
        new(r.end) value_type(+*p);
    }
    return r;
  }

  template<class T>
  Matrix<T> Matrix<T>::operator-() const {
    Matrix<value_type> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(std::malloc(sizeof(value_type) * _size));
      for (const value_type *p = begin; n; ++r.end, ++p, --n)
        new(r.end) value_type(-*p);
    }
  }

  template<class T>
  Matrix<T> Matrix<T>::operator~() const {
    Matrix<value_type> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(std::malloc(sizeof(value_type) * _size));
      for (const value_type *p = begin; n; ++r.end, ++p, --n)
        new(r.end) value_type(~*p);
    }
  }

  template<class T>
  Matrix<bool> Matrix<T>::operator!() const {
    Matrix<bool> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(std::malloc(sizeof(bool) * _size));
      for (const value_type *p = begin; n; ++r.end, ++p, --n)
        new(r.end) bool(~*p);
    }
    return r;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator*=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p *= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator/=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p /= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator%=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p %= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator+=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p += x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator-=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p -= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator^=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p ^= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator&=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p &= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator|=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p |= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator<<=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p <<= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator>>=(const value_type &x) {
    for (value_type *p = begin; p != end; ++p)
      *p >>= x;
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator*=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) *= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator/=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) /= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator%=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) %= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator+=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) += m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator-=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) -= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator^=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) ^= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator|=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) |= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator&=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) &= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator<<=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) <<= m[i];
    return *this;
  }

  template<class T>
  Matrix<T> &Matrix<T>::operator>>=(const Matrix &m) {
    for (size_t i = 0; i < m.size(); i++)
      at(i) >>= m[i];
    return *this;
  }

  template<class T>
  void Matrix<T>::clear(size_t capacity) {
    if (begin != nullptr) {
      while (end != begin)
        (--end)->~value_type();
      std::free(begin);
      begin = end = nullptr;
    }
  }

  template<class T>
  void Matrix<T>::resize(size_t n, value_type x) {
    _cols = _size = n;
    _rows = 1;
    clear(size());
    if (n) {
      begin = end = static_cast<value_type *>(std::malloc(n * sizeof(value_type)));
      try {
        for (size_t n_left = n; n_left; --n_left, ++end)
          new(end) value_type(x);
      } catch (...) {
        clear(n);
        throw;
      }
    }
  }

  template<class T>
  void Matrix<T>::resize(size_t n, size_t m, value_type x) {
    _rows = n;
    _cols = m;
    _size = n * m;
    clear(size());
    if (_size) {
      begin = end = static_cast<value_type *>(std::malloc(_size * sizeof(value_type)));
      try {
        for (size_t n_left = _size; n_left; --n_left, ++end)
          new(end) value_type(x);
      } catch (...) {
        clear(_size);
        throw;
      }
    }
  }

  template<class T>
  Matrix<T> &Matrix<T>::assign_range(const value_type *f, const value_type *l) {
    size_t n = l - f;
    if (size() != n) {
      clear(size());
      begin = static_cast<value_type *>(std::malloc(n * sizeof(value_type)));
      end = begin + n;
      std::uninitialized_copy(f, l, begin);
    } else
      std::copy(f, l, begin);
    return *this;
  }

  template<class T, class Function>
  Matrix<T> apply_operator(const Matrix<T> &x, Function fn) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = fn(x[i]);
    return output;
  }

#ifdef USES_OPENBLAS

  Matrix<float> operator*(const Matrix<float> &x, const Matrix<float> &y) {
    const int n = static_cast<int>(x.rows()), m = static_cast<int>(y.cols()), k = static_cast<int>(x.cols());
    const int lda = n, ldb = k, ldc = n;
    const float *a_ptr = x.data(), *b_ptr = y.data();
    const float alpha = 1.0, beta = 0.0;
    Matrix<float> output(x.rows(), y.cols());
    float *c_ptr = output.data();
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, k, alpha, a_ptr, lda, b_ptr, ldb, beta, c_ptr, ldc);
    return output;
  }

  Matrix<double> operator*(const Matrix<double> &x, const Matrix<double> &y) {
    const int n = static_cast<int>(x.rows()), m = static_cast<int>(y.cols()), k = static_cast<int>(x.cols());
    const int lda = n, ldb = k, ldc = n;
    const double *a_ptr = x.data(), *b_ptr = y.data();
    const double alpha = 1.0, beta = 0.0;
    Matrix<double> output(x.rows(), y.cols());
    double *c_ptr = output.data();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, k, alpha, a_ptr, lda, b_ptr, ldb, beta, c_ptr, ldc);
    return output;
  }

#endif

  template<class T>
  Matrix<T> operator*(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), y.cols());
    output = 0;
    for (size_t i = 0; i < x.rows(); i++) {
      for (size_t j = 0; j < y.cols(); j++) {
        for (size_t k = 0; k < x.cols(); k++)
          output(i, j) += x(i, k) * y(k, j);
      }
    }
    return output;
  }

#ifdef USES_OPENBLAS

  Matrix<float> operator*(const Matrix<float> &x, const float &y) {
      Matrix<float> output = x;
      cblas_sscal(static_cast<int>(output.size()), y, output.data(), 1);
      return output;
  }

  Matrix<double> operator*(const Matrix<double> &x, const double &y) {
      Matrix<double> output = x;
      cblas_dscal(static_cast<int>(output.size()), y, output.data(), 1);
      return output;
  }

#endif

  template<class T>
  Matrix<T> operator*(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] * y;
    return output;
  }

#ifdef USES_OPENBLAS

  Matrix<float> operator*(const float &x, const Matrix<float> &y) {
      return operator*(y, x);
  }

  Matrix<double> operator*(const double &x, const Matrix<double> &y) {
      return operator*(y, x);
  }

#endif

  template<class T>
  Matrix<T> operator*(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x * y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator/(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] / y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator/(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] / y;
    return output;
  }

  template<class T>
  Matrix<T> operator/(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x / y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator%(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] % y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator%(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] % y;
    return output;
  }

  template<class T>
  Matrix<T> operator%(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output = x % y[i];
    return output;
  }

#ifdef USES_OPENBLAS

  Matrix<float> operator+(const Matrix<float> &x, const Matrix<float> &y) {
    Matrix<float> output = y;
    const int n = static_cast<int>(x.size());
    const int incY = 1, incX = 1;
    const float alpha = 1.0;
    cblas_saxpy(n, alpha, x.data(), incX, output.data(), incY);
    return output;
  }

  Matrix<double> operator+(const Matrix<double> &x, const Matrix<double> &y) {
    Matrix<double> output = y;
    const int n = static_cast<int>(x.size());
    const int incY = 1, incX = 1;
    const double alpha = 1.0;
    cblas_daxpy(n, alpha, x.data(), incX, output.data(), incY);
    return output;
  }

#endif

  template<class T>
  Matrix<T> operator+(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] + y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator+(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] + y;
    return output;
  }

  template<class T>
  Matrix<T> operator+(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x + y[i];
    return output;
  }

#ifdef USES_OPENBLAS

  Matrix<float> operator-(const Matrix<float> &x, const Matrix<float> &y) {
    Matrix<float> output = x;
    const int n = static_cast<int>(x.size());
    const int incY = 1, incX = 1;
    const auto alpha = static_cast<float>(-1.0);
    cblas_saxpy(n, alpha, y.data(), incX, output.data(), incY);
    return output;
  }

  Matrix<double> operator-(const Matrix<double> &x, const Matrix<double> &y) {
    Matrix<double> output = x;
    const int n = static_cast<int>(x.size());
    const int incY = 1, incX = 1;
    const double alpha = -1.0;
    cblas_daxpy(n, alpha, y.data(), incX, output.data(), incY);
    return output;
  }

#endif

  template<class T>
  Matrix<T> operator-(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] - y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator-(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] - y;
    return output;
  }

  template<class T>
  Matrix<T> operator-(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x - y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator^(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] ^ y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator^(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] ^ y;
    return output;
  }

  template<class T>
  Matrix<T> operator^(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x ^ y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator&(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] & y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator&(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] & y;
    return output;
  }

  template<class T>
  Matrix<T> operator&(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x & y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator|(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] | y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator|(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] | y;
    return output;
  }

  template<class T>
  Matrix<T> operator|(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x | y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator<<(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] << y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator<<(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] << y;
    return output;
  }

  template<class T>
  Matrix<T> operator<<(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x << y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator>>(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] >> y[i];
    return output;
  }

  template<class T>
  Matrix<T> operator>>(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] >> y;
    return output;
  }

  template<class T>
  Matrix<T> operator>>(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x >> y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator&&(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] && y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator&&(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] && y;
    return output;
  }

  template<class T>
  Matrix<bool> operator&&(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x && y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator||(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] || y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator||(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] || y;
    return output;
  }

  template<class T>
  Matrix<bool> operator||(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x || y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator==(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] == y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator==(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] == y;
    return output;
  }

  template<class T>
  Matrix<bool> operator==(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x == y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator!=(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] != y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator!=(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] != y;
    return output;
  }

  template<class T>
  Matrix<bool> operator!=(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x != y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator<(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] < y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator<(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] < y;
    return output;
  }

  template<class T>
  Matrix<bool> operator<(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x < y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator>(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] > y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator>(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] > y;
    return output;
  }

  template<class T>
  Matrix<bool> operator>(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x > y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator<=(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] <= y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator<=(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] <= y;
    return output;
  }

  template<class T>
  Matrix<bool> operator<=(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x <= y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator>=(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] >= y[i];
    return output;
  }

  template<class T>
  Matrix<bool> operator>=(const Matrix<T> &x, const T &y) {
    Matrix<bool> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] >= y;
    return output;
  }

  template<class T>
  Matrix<bool> operator>=(const T &x, const Matrix<T> &y) {
    Matrix<bool> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = x >= y[i];
    return output;
  }

  template<class T>
  Matrix<T> abs(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::abs(in); });
  }

  template<class T>
  Matrix<T> acos(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::acos(in); });
  }

  template<class T>
  Matrix<T> asin(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::asin(in); });
  }

  template<class T>
  Matrix<T> atan2(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = std::atan2(x[i], y[i]);
    return output;
  }

  template<class T>
  Matrix<T> atan2(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = std::atan2(x[i], y);
    return output;
  }

  template<class T>
  Matrix<T> atan2(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = std::atan2(x, y[i]);
    return output;
  }

  template<class T>
  Matrix<T> cos(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::cos(in); });
  }

  template<class T>
  Matrix<T> cosh(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::cosh(in); });
  }

  template<class T>
  Matrix<T> exp(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::exp(in); });
  }

  template<class T>
  Matrix<T> log(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::log(in); });
  }

  template<class T>
  Matrix<T> log10(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::log10(in); });
  }

  template<class T>
  Matrix<T> pow(const Matrix<T> &x, const Matrix<T> &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = std::pow(x[i], y[i]);
    return output;
  }

  template<class T>
  Matrix<T> pow(const Matrix<T> &x, const T &y) {
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = std::pow(x[i], y);
    return output;
  }

  template<class T>
  Matrix<T> pow(const T &x, const Matrix<T> &y) {
    Matrix<T> output(y.rows(), y.cols());
    for (size_t i = 0; i < y.size(); i++)
      output[i] = std::pow(x, y[i]);
    return output;
  }

  template<class T>
  Matrix<T> sin(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::sin(in); });
  }

  template<class T>
  Matrix<T> sinh(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::sinh(in); });
  }

  template<class T>
  Matrix<T> sqrt(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::sqrt(in); });
  }

  template<class T>
  Matrix<T> tan(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::tan(in); });
  }

  template<class T>
  Matrix<T> tanh(const Matrix<T> &x) {
    return apply_operator(x, [](T in) { return std::tanh(in); });
  }

  template<class T>
  bool all_equal(const Matrix<T> &x, const Matrix<T> &y) {
    if (x.size() != y.size()) return false;
    if (x.rows() != y.rows()) return false;
    if (x.cols() != y.cols()) return false;
    for (size_t i = 0; i < x.size(); i++) {
      if (x[i] != y[i]) return false;
    }
    return true;
  }

  template<class T>
  Matrix<T> transpose(const Matrix<T> &x) {
    const size_t N = x.rows(), M = x.cols();
    Matrix<T> output(M, N);
    for (size_t n = 0; n < N * M; n++) {
      size_t i = n / N, j = n % N;
      output[n] = x[M * j + i];
    }
    return output;
  }

  template<class T>
  Matrix<T> kronecker_product(const Matrix<T> &x, const Matrix<T> &y) {
    size_t p = y.rows(), q = y.cols();
    size_t N = x.rows() * p, M = x.cols() * q;
    Matrix<T> output(N, M);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        output(i, j) = x(i / p, j / q) * y(i % p, j % q);
      }
    }
    return output;
  }

  template<class T>
  Matrix<T> box_product(const Matrix<T> &x, const Matrix<T> &y) {
    size_t p = y.rows(), q = y.cols();
    size_t N = x.rows() * p, M = x.cols() * q;
    Matrix<T> output(N, M);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        output(i, j) = x(i / p, j % q) * y(i % p, j / q);
      }
    }
    return output;
  }

  template<class T>
  Matrix<T> hadamard_product(const Matrix<T> &x, const Matrix<T> &y) {
    if (x.rows() != y.rows() || x.cols() != y.cols() || x.size() != y.size())
      throw;
    Matrix<T> output(x.rows(), x.cols());
    for (size_t i = 0; i < x.size(); i++)
      output[i] = x[i] * y[i];
    return output;
  }

  template<class T>
  Matrix<T> concat(CONCAT_TYPE contype, const Matrix<T> &x, const Matrix<T> &y) {
    if (contype == HORIZONTAL) {
      Matrix<T> output(x.rows(), x.cols() + y.cols());
      for (size_t i = 0; i < x.rows(); i++) {
        for (size_t j = 0; j < x.cols(); j++)
          output(i, j) = x(i, j);
        for (size_t j = x.cols(); j < x.cols() + y.cols(); j++)
          output(i, j) = y(i, j - x.cols());
      }
      return output;
    } else {
      Matrix<T> output(x.rows() + y.rows(), x.cols());
      for (size_t i = 0; i < x.rows(); i++) {
        for (size_t j = 0; j < x.cols(); j++)
          output(i, j) = x(i, j);
      }
      for (size_t i = x.rows(); i < x.rows() + y.rows(); i++) {
        for (size_t j = 0; j < x.cols(); j++)
          output(i, j) = y(i - x.rows(), j);
      }
      return output;
    }
  }
}

#endif //CPPNNET2_MATRIX_H
