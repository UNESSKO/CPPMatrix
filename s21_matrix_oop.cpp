
#include "./s21_matrix_oop.h"

S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_),
      cols_(other.cols_),
      matrix_(new double[rows_ * cols_]) {
  std::copy(other.matrix_, other.matrix_ + rows_ * cols_, matrix_);
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ < 0 || cols_ < 0) {
    throw std::length_error("Matrix size must be greater or equal than 0");
  }

  matrix_ = new double[rows_ * cols_]{};
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix copy{other};
  *this = std::move(copy);
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this != &other) {
    Free();

    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }
  return *this;
}

S21Matrix::~S21Matrix() noexcept { Free(); }

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (std::abs(other(i, j) - (*this)(i, j)) > kEpsilon) {
        return false;
      }
    }
  }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("Incorrect matrix size for Sum");
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      (*this)(i, j) += other(i, j);
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("Incorrect matrix size for Sub");
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      (*this)(i, j) -= other(i, j);
    }
  }
}

void S21Matrix::MulNumber(const double number) noexcept {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      (*this)(i, j) *= number;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.get_rows()) {
    throw std::logic_error("Incorrect matrix size for Multiplication");
  }

  S21Matrix result{rows_, other.get_cols()};

  for (int i = 0; i < result.get_rows(); ++i) {
    for (int j = 0; j < result.get_cols(); ++j) {
      for (int k = 0; k < cols_; ++k) {
        result(i, j) += (*this)(i, k) * other(k, j);
      }
    }
  }
  *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result{cols_, rows_};

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result(j, i) = (*this)(i, j);
    }
  }

  return result;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::logic_error("Incorrect matrix size for CalcComplements");
  }

  S21Matrix result{rows_, cols_};

  for (int i = 0; i < result.get_rows(); ++i) {
    for (int j = 0; j < result.get_cols(); ++j) {
      S21Matrix minor_matrix = GetMinorMatrix(i, j);
      result(i, j) = minor_matrix.Determinant();

      if ((i + j) % 2 != 0) {
        result(i, j) = -result(i, j);
      }
    }
  }

  return result;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::logic_error("Incorrect matrix size for Determinant");
  }

  double result = 1.0;
  S21Matrix tmp{*this};
  int size = rows_;

  for (int i = 0; i < size; ++i) {
    int pivoting = i;
    for (int j = i + 1; j < size; ++j) {
      if (std::abs(tmp(j, i)) > std::abs(tmp(pivoting, i))) {
        pivoting = j;
      }
    }

    if (std::abs(tmp(pivoting, i)) < kEpsilon) {
      return 0.0;
    }

    tmp.SwapRows(i, pivoting);
    result *= tmp(i, i);

    if (i != pivoting) {
      result = -result;
    }

    for (int j = i + 1; j < size; ++j) {
      double koef = tmp(j, i) / tmp(i, i);
      for (int k = i; k < size; ++k) {
        tmp(j, k) -= tmp(i, k) * koef;
      }
    }
  }

  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::logic_error("Incorrect matrix size for Inverse");
  }

  double det = Determinant();

  if (std::abs(det) < kEpsilon) {
    throw std::logic_error("Determinant must be non-zero to calculate Inverse");
  }

  return Transpose().CalcComplements() * (1.0 / det);
}

int S21Matrix::get_rows() const noexcept { return rows_; }

int S21Matrix::get_cols() const noexcept { return cols_; }

void S21Matrix::set_rows(int new_rows) {
  if (new_rows < 0) {
    throw std::length_error("matrix rows count must be non-negative");
  }

  if (new_rows != rows_) {
    S21Matrix tmp{new_rows, cols_};
    for (int i = 0; i < std::min(new_rows, rows_); ++i) {
      for (int j = 0; j < cols_; ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }
    *this = std::move(tmp);
  }
}

void S21Matrix::set_cols(int new_cols) {
  if (new_cols < 0) {
    throw std::length_error("matrix columns count must be non-negative");
  }

  if (new_cols != cols_) {
    S21Matrix tmp{rows_, new_cols};
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < std::min(new_cols, cols_); ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }
    *this = std::move(tmp);
  }
}

double &S21Matrix::operator()(int row, int col) & {
  return const_cast<double &>(GetMatrixElement(row, col));
}

const double &S21Matrix::operator()(int row, int col) const & {
  return GetMatrixElement(row, col);
}

bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix tmp{*this};
  tmp.SumMatrix(other);
  return tmp;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix tmp{*this};
  tmp.SubMatrix(other);
  return tmp;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*(double number) const noexcept {
  S21Matrix tmp{*this};
  tmp.MulNumber(number);
  return tmp;
}

S21Matrix operator*(double number, const S21Matrix &matrix) noexcept {
  S21Matrix tmp = matrix * number;
  return tmp;
}

S21Matrix &S21Matrix::operator*=(double number) {
  MulNumber(number);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix tmp{*this};
  tmp.MulMatrix(other);
  return tmp;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

void S21Matrix::Free() noexcept {
  delete[] matrix_;
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}

const double &S21Matrix::GetMatrixElement(int row, int col) const {
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Incorrect input for (), index is out of range.");
  }
  return matrix_[row * cols_ + col];
}

void S21Matrix::SwapRows(int row1, int row2) noexcept {
  if (row1 != row2) {
    for (int i = 0; i < cols_; ++i) {
      std::swap((*this)(row1, i), (*this)(row2, i));
    }
  }
}

S21Matrix S21Matrix::GetMinorMatrix(const int skip_row,
                                    const int skip_column) const {
  S21Matrix result{rows_ - 1, cols_ - 1};

  int shift_i = 0;
  for (int i = 0; i < result.get_rows(); ++i) {
    if (i == skip_row) {
      shift_i = 1;
    }
    int shift_j = 0;
    for (int j = 0; j < result.get_cols(); ++j) {
      if (j == skip_column) {
        shift_j = 1;
      }
      result(i, j) = (*this)(i + shift_i, j + shift_j);
    }
  }

  return result;
}