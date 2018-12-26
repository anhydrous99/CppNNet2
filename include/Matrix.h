//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_MATRIX_H
#define CPPNNET2_MATRIX_H

#include <cstddef>
#include <algorithm>
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

    void swap(Matrix &v) noexcept;

    constexpr size_t size();
    constexpr size_t rows();
    constexpr size_t cols();

    void clear(size_t capacity);
    void resize(size_t n, value_type x = value_type());
    void resize(size_t n, size_t m, value_type x = value_type());

  private:
    Matrix &assign_range(const value_type *f, const value_type *l);
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
    _cols = _size = n;
    _rows = 1;
    begin = p;
    end = p + n;
  }

  template<class T>
  Matrix<T>::Matrix(const Matrix<T> &m) : begin(0), end(0) {
    _rows = m.rows();
    _cols = m.cols();
    _size = m.size();
    if (m.size()) {
      begin = end = static_cast<value_type *>(malloc(sizeof(value_type) * m.size()));
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
  Matrix<T>::Matrix(size_t n, size_t m) : begin(0), end(0) {
    _rows = n;
    _cols = m;
    _size = n * m;
    if (n && m) {
      begin = end = static_cast<value_type *>(malloc(sizeof(value_type) * _size));
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
    resize(n * m, x);
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
    if (this != *m)
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
  }

  template<class T>
  Matrix<T> Matrix<T>::operator+() const {
    Matrix<value_type> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(malloc(sizeof(value_type) * _size));
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
      r.begin = r.end = static_cast<value_type *>(malloc(sizeof(value_type) * _size));
      for (const value_type *p = begin; n; ++r.end, ++p, --n)
        new(r.end) value_type(-*p);
    }
  }

  template<class T>
  Matrix<T> Matrix<T>::operator~() const {
    Matrix<value_type> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(malloc(sizeof(value_type) * _size));
      for (const value_type *p = begin; n; ++r.end, ++p, --n)
        new(r.end) value_type(~*p);
    }
  }

  template<class T>
  Matrix<bool> Matrix<T>::operator!() const {
    Matrix<bool> r;
    size_t n = size();
    if (n) {
      r.begin = r.end = static_cast<value_type *>(malloc(sizeof(bool) * _size));
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
}

#endif //CPPNNET2_MATRIX_H
