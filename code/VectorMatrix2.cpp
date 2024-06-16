#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

/*This code contains additional functionality to the vector and matrix classes*/

const char greska[] = "Incompatible formats";

class Vector {
  std::vector<double> *elementi;

public:
  explicit Vector(int n);
  Vector(std::initializer_list<double> l);
  Vector(const Vector &v) : elementi(new std::vector<double>(*v.elementi)) {}
  int NElems() const { return elementi->size(); }
  double &operator[](int i) { return (*elementi)[i]; }
  double operator[](int i) const { return (*elementi)[i]; }
  double &operator()(int i) {
    if (i <= 0 || i > NElems())
      throw std::range_error("Invalid index");
    return (*elementi)[i - 1];
  }
  double operator()(int i) const {
    if (i <= 0 || i > NElems())
      throw std::range_error("Invalid index");
    return (*elementi)[i - 1];
  }
  double Norm() const {
    double x = 0;
    for (int i = 0; i < NElems(); i++)
      x += (*elementi)[i] * (*elementi)[i];
    return sqrt(x);
  }
  friend double VectorNorm(const Vector &v) { return v.Norm(); }
  double GetEpsilon() const {
    return 10 * Norm() * std::numeric_limits<double>::epsilon();
  }
  void Print(char separator = '\n', double eps = -1) const;
  friend void PrintVector(const Vector &v, char separator = '\n',
                          double eps = -1) {
    v.Print(separator, eps);
  }
  friend Vector operator+(const Vector &v1, const Vector &v2);
  Vector &operator+=(const Vector &v);
  friend Vector operator-(const Vector &v1, const Vector &v2);
  Vector &operator-=(const Vector &v);
  friend Vector operator*(double s, const Vector &v);
  friend Vector operator*(const Vector &v, double s);
  Vector &operator*=(double s);
  friend double operator*(const Vector &v1, const Vector &v2);
  friend Vector operator/(const Vector &v, double s);
  Vector &operator/=(double s);
  ~Vector() { delete elementi; }
  Vector &operator=(const Vector &v) {
    (*elementi) = (*v.elementi);
    return *this;
  }
  void Chop(double eps = -1) {
    if (eps < 0)
      eps = GetEpsilon();
    for (int i = 0; i < NElems(); i++)
      if ((*this)[i] < eps)
        (*this)[i];
  }
  bool EqualTo(const Vector &v, double eps = -1) const {
    if (v.NElems() != NElems())
      return false;
    if (eps < 0)
      eps = GetEpsilon();
    for (int i = 0; i < NElems(); i++)
      if (fabs((*this)[i] - v[i]) >= eps)
        return false;
    return true;
  }
  Vector(Vector &&v) : elementi(v.elementi) { v.elementi = nullptr; }
};

Vector::Vector(int n) {
  if (n <= 0)
    throw std::range_error("Bad dimension");
  elementi = new std::vector<double>(n);
}

Vector::Vector(std::initializer_list<double> l) {
  if (l.size() == 0)
    throw std::range_error("Bad dimension");
  elementi = new std::vector<double>(l);
}

void Vector::Print(char separator, double eps) const {
  if (eps < 0)
    eps = GetEpsilon();
  int i;
  for (i = 0; i < NElems() - 1; i++) {
    double x = (*elementi)[i];
    if (fabs(x) < eps)
      x = 0;
    std::cout << x << separator;
  }
  std::cout << (*elementi)[i];
  if (separator == '\n')
    std::cout << '\n';
}

Vector operator+(const Vector &v1, const Vector &v2) {
  if (v1.NElems() != v2.NElems())
    throw std::domain_error(greska);
  Vector v3(v1);
  for (int i = 0; i < v1.NElems(); i++)
    (*v3.elementi)[i] += (*v2.elementi)[i];
  return v3;
}

Vector operator-(const Vector &v1, const Vector &v2) {
  if (v1.NElems() != v2.NElems())
    throw std::domain_error(greska);
  Vector v3(v1);
  for (int i = 0; i < v1.NElems(); i++)
    (*v3.elementi)[i] -= (*v2.elementi)[i];
  return v3;
}

Vector &Vector::operator+=(const Vector &v) {
  if (NElems() != v.NElems())
    throw std::domain_error(greska);
  for (int i = 0; i < NElems(); i++)
    (*elementi)[i] += (*v.elementi)[i];
  return *this;
}

Vector &Vector::operator-=(const Vector &v) {
  if (NElems() != v.NElems())
    throw std::domain_error(greska);
  for (int i = 0; i < NElems(); i++)
    (*elementi)[i] -= (*v.elementi)[i];
  return *this;
}

Vector operator*(double s, const Vector &v) {
  Vector izlaz(v);
  for (int i = 0; i < izlaz.NElems(); i++)
    (*izlaz.elementi)[i] *= s;
  return izlaz;
}

Vector operator*(const Vector &v, double s) { return s * v; }

Vector &Vector::operator*=(double s) {
  for (int i = 0; i < NElems(); i++)
    (*elementi)[i] *= s;
  return *this;
}

double operator*(const Vector &v1, const Vector &v2) {
  if (v1.NElems() != v2.NElems())
    throw std::domain_error(greska);
  double x = 0;
  for (int i = 0; i < v1.NElems(); i++)
    x += (*v1.elementi)[i] * (*v2.elementi)[i];
  return x;
}

Vector &Vector::operator/=(double s) {
  if (!s)
    throw std::domain_error("Division by zero");
  for (int i = 0; i < NElems(); i++)
    (*elementi)[i] /= s;
  return *this;
}

Vector operator/(const Vector &v, double s) {
  Vector k(v);
  k /= s;
  return k;
}

class Matrix {
  std::vector<std::vector<double>> elemente;

public:
  friend Matrix LeftDiv(Matrix m1, Matrix m2);
  friend Vector LeftDiv(Matrix m, Vector v);
  friend Matrix operator/(Matrix m, double s);
  Matrix &operator/=(double s);
  friend Matrix operator/(Matrix m1, Matrix m2);
  Matrix &operator/=(Matrix m);
  double Det() const;
  friend double Det(Matrix m);
  void Invert();
  friend Matrix Inverse(Matrix m);
  void ReduceToRREF();
  friend Matrix RREF(Matrix m);
  int Rank() const;
  friend int Rank(Matrix m);
  void Chop(double eps = -1) {
    if (eps < 0)
      eps = GetEpsilon();
    for (int i = 0; i < NRows(); i++) {
      for (int j = 0; j < NCols(); j++) {
        (*this)[i][j] = 0;
      }
    }
  }

  bool EqualTo(const Matrix &m, double eps = -1) const {
    if (eps < 0)
      eps = GetEpsilon();
    if (NRows() != m.NRows() || NCols() != m.NCols())
      return false;
    for (int i = 0; i < NRows(); i++) {
      for (int j = 0; j < NCols(); j++)
        if (fabs((*this)[i][j] - m[i][j]) >= eps)
          return false;
    }
    return true;
  }

  Matrix(int m, int n) {
    if (m <= 0 || n <= 0)
      throw std::range_error("Bad dimension");
    elemente = std::vector<std::vector<double>>(m, std::vector<double>(n));
  }
  Matrix(const Vector &v) {
    elemente =
        std::vector<std::vector<double>>(v.NElems(), std::vector<double>(1));
    for (int i = 0; i < v.NElems(); i++)
      elemente[i][0] = v[i];
  }
  Matrix(std::initializer_list<std::vector<double>> l);
  int NRows() const { return elemente.size(); }
  int NCols() const {
    if (!elemente.size())
      return 0;
    return elemente.at(0).size();
  }
  double *operator[](int i) { return elemente[i].data(); }
  const double *operator[](int i) const { return elemente[i].data(); }
  double &operator()(int i, int j) {
    if (i <= 0 || i > NRows() || j <= 0 || j > NCols())
      throw std::range_error("Invalid index");
    return elemente[i - 1][j - 1];
  }
  double operator()(int i, int j) const {
    if (i <= 0 || i > NRows() || j <= 0 || j > NCols())
      throw std::range_error("Invalid index");
    return elemente[i - 1][j - 1];
  }
  double Norm() const;
  friend double MatrixNorm(const Matrix &m) { return m.Norm(); }
  double GetEpsilon() const {
    return 10 * Norm() * std::numeric_limits<double>::epsilon();
  }
  void Print(int width = 10, double eps = -1) const;
  friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {
    m.Print(width, eps);
  }
  /*nesto nesto opasno friend definisat u klasu*/
  friend Matrix operator+(Matrix m1, const Matrix &m2) { return m1 += m2; }
  Matrix &operator+=(const Matrix &m);
  friend Matrix operator-(Matrix m1, const Matrix &m2) { return m1 -= m2; }
  Matrix &operator-=(const Matrix &m);
  friend Matrix operator*(double s, const Matrix &m) { return m * s; }
  friend Matrix operator*(Matrix m, double s) { return m *= s; }
  Matrix &operator*=(double s);
  friend Matrix operator*(Matrix m1, const Matrix &m2) { return m1 *= m2; }
  Matrix &operator*=(const Matrix &m);
  friend Vector operator*(const Matrix &m, const Vector &v);
  friend Matrix Transpose(Matrix m) {
    m.Transpose();
    return m;
  }
  void Transpose();
};

double Matrix::Norm() const {
  double x = 0;
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      x += (*this)[i][j] * (*this)[i][j];
    }
  }
  return sqrt(x);
}

void Matrix::Print(int width, double eps) const {
  if (eps < 0)
    eps = GetEpsilon();
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      double x = (*this)[i][j];
      if (fabs(x) < eps)
        x = 0;
      std::cout << std::right << std::setw(width) << x;
    }
    std::cout << '\n';
  }
}

Matrix &Matrix::operator+=(const Matrix &m) {
  if (NRows() != m.NRows() || NCols() != m.NCols())
    throw std::domain_error(greska);
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      (*this)[i][j] += m[i][j];
    }
  }
  return *this;
}

Matrix &Matrix::operator-=(const Matrix &m) {
  if (NRows() != m.NRows() || NCols() != m.NCols())
    throw std::domain_error(greska);
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      (*this)[i][j] -= m[i][j];
    }
  }
  return *this;
}

Matrix &Matrix::operator*=(double s) {
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      (*this)[i][j] *= s;
    }
  }
  return *this;
}

Matrix &Matrix::operator*=(const Matrix &m) {
  if (NCols() != m.NRows())
    throw std::domain_error(greska);
  Matrix help = (*this);
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      double c = 0;
      for (int k = 0; k < NCols(); k++)
        c += (*this)[i][k] * m[k][j];
      help[i][j] = c;
    }
  }
  elemente = std::move(help.elemente);
  return *this;
}

Vector operator*(const Matrix &m, const Vector &v) {
  if (v.NElems() != m.NCols())
    throw std::domain_error(greska);
  Vector pomocni(m.NRows());
  for (int i = 0; i < pomocni.NElems(); i++)
    for (int j = 0; j < m.NCols(); j++)
      pomocni(i + 1) += m(i + 1, j + 1) * v(j + 1);
  return pomocni;
}

void Matrix::Transpose() {
  if (NRows() == NCols())
    for (int i = 0; i < NRows(); i++)
      for (int j = i + 1; j < NCols(); j++)
        std::swap((*this)[i][j], (*this)[j][i]);
  else {
    std::vector<std::vector<double>> pomocna(NCols(),
                                             std::vector<double>(NRows()));
    for (int i = 0; i < NRows(); i++) {
      for (int j = 0; j < NCols(); j++) {
        pomocna[j][i] = (*this)[i][j];
      }
    }
    elemente = std::move(pomocna);
  }
}

Matrix::Matrix(std::initializer_list<std::vector<double>> l) {
  if (!l.size() || !l.begin()->size())
    throw std::range_error("Bad dimension");
  int k = (*l.begin()).size(), i = 0;
  std::vector<std::vector<double>> pomoc(l.size(), std::vector<double>(k));
  for (const std::vector<double> *p = l.begin(); p != l.end(); p++, i++) {
    if (p->size() != k)
      throw std::logic_error("Bad matrix");
    pomoc[i] = *p;
  }
  elemente = std::move(pomoc);
}

Matrix LeftDiv(Matrix m1, Matrix m2) {
  int n = m1.NRows(), m = m2.NCols();
  if (m1.NCols() != n)
    throw std::domain_error("Divisor matrix is not square");
  if (n != m2.NRows())
    throw std::domain_error("Incompatible formats");
  for (int k = 0; k < n; k++) {
    int p = k;
    for (int i = k + 1; i < n; i++)
      if (fabs(m1[i][k]) > fabs(m1[p][k]))
        p = i;
    if (fabs(m1[p][k]) < m1.GetEpsilon())
      throw std::domain_error("Divisor matrix is singular");
    if (p != k) {
      m1.elemente[p].swap(m1.elemente[k]);
      m2.elemente[p].swap(m2.elemente[k]);
    }
    for (int i = k + 1; i < n; i++) {
      double mi = m1[i][k] / m1[k][k];
      for (int j = k + 1; j < n; j++)
        m1[i][j] -= mi * m1[k][j];
      for (int j = 0; j < m; j++)
        m2[i][j] -= mi * m2[k][j];
    }
  }
  for (int k = 0; k < m; k++) {
    for (int i = m2.NRows() - 1; i > -1; i--) {
      double s = m2[i][k];
      for (int j = i + 1; j < m2.NRows(); j++)
        s -= m1[i][j] * m2[j][k];
      m2[i][k] = s / m1[i][i];
    }
  }
  return m2;
}

Vector LeftDiv(Matrix m, Vector v) {
  int n = m.NRows(), b = v.NElems();
  if (m.NCols() != n)
    throw std::domain_error("Divisor matrix is not square");
  if (n != b)
    throw std::domain_error("Incompatible formats");
  for (int k = 0; k < n; k++) {
    int p = k;
    for (int i = k + 1; i < n; i++)
      if (fabs(m[i][k]) > fabs(m[p][k]))
        p = i;
    if (fabs(m[p][k]) < m.GetEpsilon())
      throw std::domain_error("Divisor matrix is singular");
    if (p != k) {
      m.elemente[p].swap(m.elemente[k]);
      std::swap(v[p], v[k]);
    }
    for (int i = k + 1; i < n; i++) {
      double mi = m[i][k] / m[k][k];
      for (int j = k + 1; j < n; j++)
        m[i][j] -= mi * m[k][j];
      v[i] -= mi * v[k];
    }
  }
  for (int i = b - 1; i > -1; i--) {
    double s = v[i];
    for (int j = i + 1; j < b; j++)
      s -= m[i][j] * v[j];
    v[i] = s / m[i][i];
  }
  return v;
}

Matrix &Matrix::operator/=(double s) {
  if (!s)
    throw std::domain_error("Division by zero");
  int n = NRows(), m = NCols();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++)
      (*this)[i][j] /= s;
  }
  return *this;
}

Matrix operator/(Matrix m, double s) { return m /= s; }

Matrix operator/(Matrix m1, Matrix m) {
    std::swap(m1,m);
  int n = m1.NCols(), b = m.NRows();
  if (m.NRows() != b)
    throw std::domain_error("Divisor matrix is not square");
  if (n != m.NCols())
    throw std::domain_error("Incompatible formats");
  for (int k = 0; k < m1.NRows(); k++) {
    int p = k;
    for (int i = k + 1; i < n; i++)
      if (fabs(m1[k][i]) > fabs(m1[k][p]))
        p = i;
    if (fabs(m1[k][p]) < m1.GetEpsilon())
      throw std::domain_error("Divisor matrix is singular");
    if (p != k) {
      for (int t = 0; t < m1.NRows(); t++) {
        std::swap(m1[t][p], m1[t][k]);
      }
      for (int t = 0; t < b; t++)
        std::swap(m[t][p], m[t][k]);
    }
    for (int i = k + 1; i < n; i++) {
      double mi = m1[k][i] / m1[k][k];
      for (int j = k + 1; j < m1.NRows(); j++)
        m1[j][i] -= mi * m1[j][k];
      for (int j = 0; j < b; j++)
        m[j][i] -= mi * m[j][k];
    }
  }
  for (int k = 0; k < b; k++) {
    for (int i = m1.NCols() - 1; i > -1; i--) {
      double s = m[k][i];
      for (int j = i + 1; j < m.NCols(); j++)
        s -= m1[j][i] * m[k][j];
      m[k][i] = s / m1[i][i];
    }
  }
  return m;
}

Matrix &Matrix::operator/=(Matrix m) {
  return (*this) = std::move((*this) / m);
}

double Matrix::Det() const {
  double d = 1;
  int i, j, k, p, n = NRows(), m = NCols();
  if (n != m)
    throw std::domain_error("Matrix is not square");
  Matrix mat = (*this);
  for (k = 0; k < n; k++) {
    p = k;
    for (int i = k + 1; i < n; i++)
      if (fabs(mat[i][k]) > fabs(mat[p][k]))
        p = i;
    if (fabs(mat[p][k]) < mat.GetEpsilon())
      return 0;
    if (p != k) {
      mat.elemente[k].swap(mat.elemente[p]);
      d *= -1;
    }
    d *= mat[k][k];
    for (i = k + 1; i < n; i++) {
      double mi = mat[i][k] / mat[k][k];
      for (j = k + 1; j < n; j++)
        mat[i][j] -= mi * mat[k][j];
    }
  }
  return d;
}

double Det(Matrix m) { return m.Det(); }

void Matrix::Invert() {
  int i, j, k, n = NRows(), m = NCols(), p;
  std::vector<int> w(n);
  if (m != n)
    throw std::domain_error("Matrix is not square");
  Matrix &a = (*this);
  for (k = 0; k < n; k++) {
    p = k;
    for (i = k + 1; i < n; i++)
      if (fabs(a[i][k]) > fabs(a[p][k]))
        p = i;
    if (fabs(a[p][k]) < a.GetEpsilon())
      throw std::domain_error("Matrix is singular");
    if (p != k)
      a.elemente[p].swap(a.elemente[k]);
    w[k] = p;
    double mi = a[k][k];
    a[k][k] = 1;
    for (j = 0; j < n; j++)
      a[k][j] /= mi;
    for (i = 0; i < n; i++) {
      if (i == k)
        continue;
      mi = a[i][k];
      a[i][k] = 0;
      for (j = 0; j < n; j++)
        a[i][j] -= mi * a[k][j];
    }
  }
  for (j = n - 1; j > -1; j--) {
    p = w[j];
    if (p == j)
      continue;
    for (i = 0; i < n; i++)
      std::swap(a[i][p], a[i][j]);
  }
}

Matrix Inverse(Matrix m) {
  m.Invert();
  return m;
}

void Matrix::ReduceToRREF() {
  int k = -1, l = -1, j, p, i, m = NCols(), n = NRows();
  double v;
  Matrix &a = (*this);
  std::vector<bool> w(n, false);
  while (k < m && l < n) {
    l++;
    k++;
    v = 0;
    while (v < GetEpsilon() && l < n) {
      p = k;
      for (i = k; i < m; i++) {
        if (fabs(a[i][l]) > v) {
          v = fabs(a[i][l]);
          p = i;
        }
      }
      if (v < GetEpsilon())
        l++;
    }
  }
  if (l < n) {
    w[l] = true;
    if (p != k)
      a.elemente[k].swap(a.elemente[p]);
    double mi = a[k][l];
    for (j = l; j < n; j++)
      a[k][j] /= mi;
    for (i = 0; i < m; i++) {
      if (i == k)
        continue;
      mi = a[i][l];
      for (int j = l; j < n; j++)
        a[i][j] -= mi * a[k][j];
    }
  }
}

Matrix RREF(Matrix m) {
  m.ReduceToRREF();
  return m;
}

int Matrix::Rank() const {
  int k = -1, l = -1, j, p, i, m = NCols(), n = NRows();
  double v;
  Matrix a = (*this);
  std::vector<bool> w(n, false);
  while (k < m && l < n) {
    l++;
    k++;
    v = 0;
    while (v < GetEpsilon() && l < n) {
      p = k;
      for (i = k; i < m; i++) {
        if (fabs(a[i][l]) > v) {
          v = fabs(a[i][l]);
          p = i;
        }
      }
      if (v < GetEpsilon())
        l++;
    }
  }
  if (l < n) {
    w[l] = true;
    if (p != k)
      a.elemente[k].swap(a.elemente[p]);
    double mi = a[k][l];
    for (j = l; j < n; j++)
      a[k][j] /= mi;
    for (i = 0; i < m; i++) {
      if (i == k)
        continue;
      mi = a[i][l];
      for (int j = l; j < n; j++)
        a[i][j] -= mi * a[k][j];
    }
  }
  return k + 1;
}

int Rank(Matrix a) {
  int k = -1, l = -1, j, p = 0, i, m = a.NCols(), n = a.NRows();
  double v;
  std::vector<bool> w(n, false);
  while (k < m && l < n) {
    l++;
    k++;
    v = 0;
    while (v < a.GetEpsilon() && l < n) {
      p = k;
      for (i = k; i < m; i++) {
        if (fabs(a[i][l]) > v) {
          v = fabs(a[i][l]);
          p = i;
        }
      }
      if (v < a.GetEpsilon())
        l++;
    }
  }
  if (l < n) {
    w[l] = true;
    if (p != k)
      a.elemente[k].swap(a.elemente[p]);
    double mi = a[k][l];
    for (j = l; j < n; j++)
      a[k][j] /= mi;
    for (i = 0; i < m; i++) {
      if (i == k)
        continue;
      mi = a[i][l];
      for (int j = l; j < n; j++)
        a[i][j] -= mi * a[k][j];
    }
  }
  return k + 1;
}

/*
Matrix &Matrix::operator /=(Matrix m) {
    int n = NCols(), b = m.NRows();
  if (m.NRows() != n)
    throw std::domain_error("Divisor matrix is not square");
  if (n != m.NCols())
    throw std::domain_error("Incompatible formats");
  for (int k = 0; k < n; k++) {
    int p = k;
    for (int i = k + 1; i < n; i++)
      if (fabs((*this)[k][i]) > fabs((*this)[k][p]))
        p = i;
    if (fabs((*this)[k][p]) < (*this).GetEpsilon())
      throw std::domain_error("Divisor matrix is singular");
    if (p != k) {
      (*this).elemente[p].swap((*this).elemente[k]);
      m.elemente[p].swap(m.elemente[k]);
    }
    for (int i = k + 1; i < n; i++) {
      double mi = (*this)[k][i] / (*this)[k][k];
      for (int j = k + 1; j < n; j++)
        (*this)[j][i] -= mi * (*this)[j][k];
      for (int j = 0; j < b; j++)
        m[j][i] -= mi * m[j][k];
    }
  }
  for (int k = 0; k < b; k++) {
    for (int i = m.NCols() - 1; i > -1; i--) {
      double s = m[k][i];
      for (int j = i + 1; j < m.NCols(); j++)
        s -= (*this)[j][i] * m[k][j];
      m[k][i] = s / (*this)[i][i];
    }
  }
  return m;
}*/

int main() {
  Matrix A{{0, 3, 2}, {4, 6, 1}, {3, 1, 7}};
  Matrix b{{4, 1, 5}, {1, 2, 1}};
  Matrix rez = b * Inverse(A);
  b / A;
  std::cout << rez.EqualTo(b / A) << std::endl;
  Matrix C{{2, 3, 5}, {2, 3, 7}, {4, 1, 8}};
  rez = b * Inverse(C);
  std::cout << rez.EqualTo(b / C);
  return 0;
}