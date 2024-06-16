#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <utility>
#include <vector>


/* A recursive helper function for the Adaptive Simpson implementation*/


template <typename FunType>
std::pair<double, bool> AdaptiveAux(FunType f, double a, double b, double eps,
                                    int maxdepth, double f1, double f2,
                                    double f3) {
  double c = (a + b) / 2;
  double I1 = (b - a) * (f1 + 4 * f3 + f2) / 6;

  double f4, f4_rjesenje = f((a + c) / 2);
  if (f4_rjesenje == std::numeric_limits<double>::infinity() ||
      f4_rjesenje == -std::numeric_limits<double>::infinity()) {
    f4 = 0;
  } else {
    f4 = f4_rjesenje;
  }

  double f5, f5_rjesenje = f((c + b) / 2);
  if (f5_rjesenje == std::numeric_limits<double>::infinity() ||
      f5_rjesenje == -std::numeric_limits<double>::infinity()) {
    f5 = 0;
  } else {
    f5 = f5_rjesenje;
  }

  double I2 = (b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12;

  if (std::abs(I1 - I2) <= eps)
    return std::pair<double, bool>(I2, true);

  if (maxdepth <= 0)
    return std::pair<double, bool>(I2, false);

  std::pair<double, bool> rjesenje1 =
      AdaptiveAux(f, a, c, eps, maxdepth - 1, f1, f3, f4);
  std::pair<double, bool> rjesenje2 =
      AdaptiveAux(f, c, b, eps, maxdepth - 1, f3, f2, f5);

  if (rjesenje1.second == false || rjesenje2.second == false) {
    return std::pair<double, bool>(rjesenje1.first + rjesenje2.first, false);
  } else {
    return std::pair<double, bool>(rjesenje1.first + rjesenje2.first, true);
  }
}


 
template <typename FunType>
std::pair<double, bool> AdaptiveIntegration(FunType f, double a, double b,
                                            double eps = 1e-10,
                                            int maxdepth = 30, int nmin = 1) {
  if (eps <= 0 || nmin <= 0 || maxdepth <= 0) {
    throw std::domain_error("Bad parameter");
  }
  int interval = 1;
  if (a > b) {
    std::swap(a, b);
    interval = -1;
  }

  std::pair<double, bool> s(0, true);
  double h = (b - a) / nmin;
  for (int i = 1; i <= nmin; i++) {
    double fa, fa_rjesenje = f(a);
    double fah, fah_rjesenje = f(a + h);
    double fah2, fah2_rjesenje = f(a + h / 2);

    if (fa_rjesenje == std::numeric_limits<double>::infinity() ||
        fa_rjesenje == -std::numeric_limits<double>::infinity()) {
      fa = 0;
    } else {
      fa = fa_rjesenje;
    }

    if (fah_rjesenje == std::numeric_limits<double>::infinity() ||
        fah_rjesenje == -std::numeric_limits<double>::infinity()) {
      fah = 0;
    } else {
      fah = fah_rjesenje;
    }

    if (fah2_rjesenje == std::numeric_limits<double>::infinity() ||
        fah2_rjesenje == -std::numeric_limits<double>::infinity()) {
      fah2 = 0;
    } else {
      fah2 = fah2_rjesenje;
    }

    std::pair<double, bool> rjesenje =
        AdaptiveAux(f, a, a + h, eps, maxdepth, fa, fah, fah2);
    s.first += rjesenje.first;

    if (rjesenje.second == false) {
      s.second = false;
    }
    a += h;
  }

  return std::pair<double, bool>(s.first * interval, s.second);
}


int main() { }