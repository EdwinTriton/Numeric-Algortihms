#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>



template <typename FunType>
std::pair<double, bool>
TanhSinhIntegration(FunType f, double a, double b, double eps = 1e-8,
                    int nmax = 1000000, int nmin = 20, double range = 3.5) {
  int interval = 1;
  if (nmax < nmin || nmax < 0 || nmin < 0 || eps < 0 || range <= 0) {
    throw std::domain_error("Bad parameter");
  }

  if (a > b) {
    std::swap(a, b);
    interval = -1;
  }

  std::function<double(double)> smjena = [a, b, f](double t) {
    double pomocna =
        f((b + a) / 2 + (b - a) / 2 * tanh(2 * std::atan(1) * sinh(t)));
    if (isnan(pomocna)) {
      pomocna = 1;
    }

    double rezultat = (((b - a) * 4 * std::atan(1) * cosh(t)) /
                       (4 * std::pow(cosh(2 * std::atan(1) * sinh(t)), 2))) *
                      pomocna;
    return rezultat;
  };

  double rb = range;
  double ra = -range;

  double N = 2;
  double h = (rb - ra) / N;
  double s, rb1, ra1;

  if (smjena(rb) == std::numeric_limits<double>::infinity() ||
      smjena(rb) == -std::numeric_limits<double>::infinity()) {
    rb1 = 0;
  } else {
    rb1 = smjena(rb);
  }

  if (smjena(ra) == std::numeric_limits<double>::infinity() ||
      smjena(ra) == -std::numeric_limits<double>::infinity()) {
    ra1 = 0;
  } else {
    ra1 = smjena(ra);
  }

  s = (ra1 + rb1) / 2;

  double Iold = s;
  double I;

  while (N < nmax) {
    for (int j = 1; j <= N / 2; j++) {
      double rjesenje = smjena(ra + (2 * j - 1) * h);
      if (rjesenje != std::numeric_limits<double>::infinity() &&
          rjesenje != -std::numeric_limits<double>::infinity()) {
        s += rjesenje;
      }
    }

    I = h * s;
    if (std::abs(I - Iold) <= eps && N >= nmin) {
      return std::pair<double, bool>(interval * I, true);
    }

    Iold = I;
    h /= 2;
    N *= 2;
  }
  return std::pair<double, bool>(interval * I, false);
}

int main() { }
