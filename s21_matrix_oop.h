//
// Created by KWAZAR_ on 28.06.2024.
//

#ifndef CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
#define CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

class S21Matrix final {
 public:
  S21Matrix() noexcept;
  explicit S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  ~S21Matrix() noexcept;

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double number) noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  int get_cols() const noexcept;
  int get_rows() const noexcept;
  void set_rows(int new_rows);
  void set_cols(int new_cols);

  double& operator()(int row, int col) &;
  double& operator()(int row, int col) && = delete;
  const double& operator()(int row, int col) const&;
  const double& operator()(int row, int col) const&& = delete;
  bool operator==(const S21Matrix& other) const;
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator*(double number) const noexcept;
  friend S21Matrix operator*(double number, const S21Matrix& matrix) noexcept;
  S21Matrix& operator*=(double number);
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix& operator*=(const S21Matrix& other);

 private:
  void Free() noexcept;
  const double& GetMatrixElement(int row, int col) const;
  void SwapRows(int row1, int row2) noexcept;
  S21Matrix GetMinorMatrix(const int skip_row, const int skip_column) const;

  int rows_;
  int cols_;
  double* matrix_;
  const double kEpsilon = 1e-7;
};

#endif  // CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
