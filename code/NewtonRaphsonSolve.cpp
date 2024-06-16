#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <functional>
#include <vector>
#include <iostream>
#include <complex>

/*Contains implementation of Regula Falsi and it's different implementations*/
enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template<typename funkcija>
bool fun(funkcija f, double x0, double &a, double &b, double hinit,double hmax, double lambda) {
    a = x0;
    double fa = f(a), h = hinit, f1 = f(a), f2;
    while (std::fabs(h) < hmax) {
        b = a+h;
        f2 = f(b);
        while (std::isnan(f2)) {
            h = h/(2*(1+lambda));
            if (fabs(h) <= fabs(a)*	std::numeric_limits<double>::epsilon()) return false;
            b = a+h;
            f2 = f(b);
        }
        if (f1*f2<=0) {
            if (b < a) std::swap(a,b);
            return true; 
        }
        h = lambda*h;
        a = b;
        f1 = f2;
    }
    return false;
}

template<typename broj>
int sgn(const broj &a) {
return (a>0) - (a<0);    
}

std::complex<double> RandomComplex() {
    double a = rand(), b = rand();
    if (int(a)%2) a*=-1;
    if (int(b)%2) a*=-1;
    while (std::fabs(a) > 10) a /=10;
    while (std::fabs(b) > 10) b /=10;
    return {a,b};
}

template <typename FunType>
 bool BracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5,double hmax = 1e10, double lambda = 1.4) {
    if (hinit < 0 || hmax < 0 || lambda < 0) throw std::domain_error("Invalid parameters"); 
    if (fun(f,x0,a,b,hinit,hmax,lambda)) return true;
    return fun(f,x0,a,b,-hinit,hmax,lambda);
 }

template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b,RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100) {
    if (eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    if (f(a)*f(b) > 0) throw std::range_error("Root must be bracketed");
    if (a > b) std::swap(a,b);
    std::function<double(double)> fi = [&] (double x) {return f(x);};
    if (mode > Illinois) {
        fi = [&] (double x) {return f(x)/(1+std::fabs(f(x)));};
    }
    double f1 = fi(a), f2 = fi(b), f3, c = a, d = b;
    int i = 0;
    while (std::fabs(c-d) > eps && i < maxiter) {
        d = c;
        c = (a*f2-b*f1)/(f2-f1);
        f3 = fi(c); 
        if (!f3) return c;
        if (f1*f3 < 0) {
            b = a;
            f2 = f1;
        }
        else if (mode%2) f2 /= 2;
        a = c;
        f1 = f3;
        i++;
    }
    if (i == maxiter) throw std::domain_error("Given accuracy has not achieved");
    return c;
}

template <typename FunType>
double RiddersSolve(FunType f, double a, double b, double eps = 1e-10,int maxiter = 100) {
    if (eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    double f1 = f(a), f2 = f(b), f3, f4, c, d;
    if (f1*f2>0) throw std::range_error("Root must be bracketed");
    int i = 0;
    while (std::fabs(b-a) > eps && i<maxiter) {
        c = (a+b)/2;
        f3 = f(c);
        if (!f3) return c;
        d = c+f3*(c-a)*sgn(f1-f2)/std::sqrt(f3*f3-f1*f2);
        f4 = f(d);
        if (!f4) return d;
        if (f3*f4 < 0) {
            a = c;
            b = d;
            f1 = f3;
            f2 = f4;
        }
        else if (f1*f4 <= 0) {
            b = d;
            f2 = f4;
        }
        else {
            a = d;
            f1 = f4;
        }
        i++;
    }
    if (i == maxiter) throw std::domain_error("Given accuracy has not achieved");
    return (a+b)/2;
}


template <typename FunType1, typename FunType2>
double NewtonRaphsonSolve(FunType1 f, FunType2 fprim, double x0, double eps = 1e-10, double damping = 0, int maxiter = 100) {
    if (eps < 0 || maxiter < 0 || damping < 0 || damping >= 1) throw std::domain_error("Invalid parameters");
    if (!damping) damping=0.5;
    double delta = eps+1, v = f(x0), d = fprim(x0), w;
    int i = 0;
    while (std::fabs(delta) > eps && i<maxiter) {
        if (std::fabs(v) <= eps) return x0;
        delta = v/d;
        w = v;
        v = f(x0-delta);
        d = fprim(x0-delta);
        while ((std::fabs(v) > std::fabs(w) || std::isnan(v) || !d) && i++ < maxiter) {
            delta = damping * delta;
            v = f(x0-delta);
            d = f(x0-delta);
        }
        x0-=delta;
    }
    if (i >= maxiter) throw std::logic_error("Convergence has not achieved");
    return x0;
}

std::pair<std::complex<double>, bool> Lag(std::vector<std::complex<double>> p, int n, std::complex<double> x, double eps, int max) {
    int k = 0;
    std::complex<double> delta = eps+1;
    while (std::fabs(delta) > eps && k < max) {
        std::complex<double> f = p.at(k), d = 0, s = 0, r;
        for (int i = n-1; i>=0; i--) {
            s=s*x+std::complex<double>(2)*d;
            d=d*x+f;
            f=f*x+p.at(i);
        }
        if (std::fabs(f) < eps) return {x,true};
        r = std::sqrt(std::complex<double>(n-1)*((std::complex<double>)(n-1)*d*d-(std::complex<double>)n*f+s));
        if (std::fabs(d+r)>std::fabs(d-r)) delta = std::complex<double>(n)*f/(d+r);
        else delta = std::complex<double>(n)*f/(d-r);
        x -= delta;
        k++;
        if (std::fabs(delta) <= eps) return {x,true};
    }
    return {x,false};
}

std::pair<std::complex<double>, bool> Lag(std::vector<double> p, int n, std::complex<double> x, double eps, int max) {
    int k = 0;
    std::complex<double> delta = eps+1;
    while (std::fabs(delta) > eps && k < max) {
        std::complex<double> f = p.at(k), d = 0, s = 0, r;
        for (int i = n-1; i>=0; i--) {
            s=s*x+std::complex<double>(2)*d;
            d=d*x+f;
            f=f*x+p.at(i);
        }
        if (std::fabs(f) < eps) return {x,true};
        r = std::sqrt(std::complex<double>(n-1)*((std::complex<double>)(n-1)*d*d-(std::complex<double>)n*f+s));
        if (std::fabs(d+r)>std::fabs(d-r)) delta = std::complex<double>(n)*f/(d+r);
        else delta = std::complex<double>(n)*f/(d-r);
        x -= delta;
        k++;
        if (std::fabs(delta) <= eps) return {x,true};
    }
    return {x,false};
}

std::vector<std::complex<double>> PolyRoots(std::vector<double> coefficients, double eps = 1e-10,int maxiters = 100, int maxtrials = 10) {
    if (eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    int i = coefficients.size();
    std::complex<double> x;
    std::vector<std::complex<double>> rijes(i);
    double v,w,u;
    i--;
    while (i>=0) { /*pazi*/
        int t = 1;
        bool c = false;
        while (!c && t<maxtrials) {
            x = RandomComplex();
            std::tie(x,c) = Lag(coefficients,i,x,eps,maxiters);
            t++;
        }
        if (!c) throw std::logic_error("Convergence has not achieved");
        if (std::imag(x)<=eps) {
            rijes.at(i) = x;
            v = coefficients.at(i);
            for (int j = i-1; j>=0; j--) {
                w = coefficients.at(j);
                coefficients.at(j) = v;
                v = w + std::real(x)*v;
            }
            i--;
        }
        else {
            rijes.at(i)=x;
            rijes.at(i-1)=std::conj(x);
            double alfa = 2*std::real(x), beta = std::fabs(x)*std::fabs(x);
            u = coefficients.at(i);
            v = coefficients.at(i-1) + alfa*u;
            for (int j = i-2; j>=0; j--) {
                w=coefficients.at(j);
                coefficients.at(j)=u;
                u=v;
                v = w+alfa*v-beta*coefficients.at(j);
            }
            i-=2;
        }
    }
    return rijes;
}

int main() {
    return 0;
}