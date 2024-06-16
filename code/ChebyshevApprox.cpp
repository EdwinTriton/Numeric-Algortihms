#include <iostream>
#include <cmath>
#include <vector>

const double pi = std::atan(1) * 4;

class ChebyshevApproximation {
    std::vector<double> c; 
    ChebyshevApproximation(std::vector<double> r, double a, double b) {m=r.size() - 1; c=std::move(r); min=a;max=b;}
    double min,max;
    int m;
    double DajMin() const {return min;}
    double DajMax() const {return max;}
    public:
    template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n);
    double operator()(double x) const;
    void set_m(int m);
    void trunc(double eps);
    double derivative(double x) const;
    ChebyshevApproximation derivative() const;
    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

template <typename FunType> 
ChebyshevApproximation::ChebyshevApproximation(FunType f, double xmin, double xmax, int n):c(n+1) {
    if (xmin>=xmax || n<1) throw std::domain_error("Bad parameters");
    min=xmin;
    max=xmax;
    m=n;
    std::vector<double> v(n+1),w(n+2);
    for (int i=0; i<=n+1; i++) w.at(i) = std::cos(pi*i/(2*n+2));
    for (int i=0; i<=std::trunc(n/2); i++) v.at(i) = f((xmax+xmin+(xmax-xmin)*w.at(2*i+1))/2);
    for (int i=std::trunc(n/2)+1; i<=n; i++) v.at(i) = f((xmin+xmax-(xmax-xmin)*w.at(2*n+1-2*i))/2);
    for (int k = 0; k <= n; k++) {
        double s = 0; 
        for (int i=0; i<=n; i++) {
            int p = (k*(2*i+1))%(4*n+4);
            if (p>2*n+2) p=4*n+4-p;
            if (p>n+1) s-=v.at(i)*w.at(2*n+2-p);
            else s+=v.at(i)*w.at(p);
        }
        c.at(k)=2*s/(n+1);
    }
}


double ChebyshevApproximation::operator()(double x) const {
    double xmin = DajMin(), xmax = DajMax(), t = (2*x-xmin-xmax)/(xmax-xmin);
    if (x<xmin || x>xmax) throw std::domain_error("Bad argument");
    double p=1,q=t,s=c.at(0)/2+c.at(1)*t;
    for (int k = 2; k<=m; k++) {
        double r = 2*t*q-p;
        s+=c.at(k)*r;
        p=q;
        q=r;
    }
    return s;
}

void ChebyshevApproximation::set_m(int m) {
    if (m<=1 || m>ChebyshevApproximation::m) throw std::domain_error("Bad order");
    ChebyshevApproximation::m = m;
}

void ChebyshevApproximation::trunc(double eps) {
    if (eps<0) throw std::domain_error("Bad tolerance");
    int i;
    for (i = m; i>=0; i--) if (std::fabs(c.at(i)) > eps) break;
    if (i<0) throw std::domain_error("Bad tolerance");
    m = i;
}

ChebyshevApproximation ChebyshevApproximation::derivative() const {
    double mi = 4/(DajMax()-DajMin());
    ChebyshevApproximation izvod = *this;
    izvod.c.at(m-1) = mi*m*c.at(m);
    izvod.c.at(m-2) = mi*(m-1)*c.at(m-1);
    for (int k = m-3; k>=0; k--) izvod.c.at(k) = izvod.c.at(k+2) + mi*(k+1)*c.at(k+1);
    izvod.m--;
    return izvod;
}

ChebyshevApproximation ChebyshevApproximation::antiderivative() const {
    std::vector<double> covi(m+2);
    int i;
    covi.at(0) = 0;
    for (i = 1; i<m; i++) covi.at(i) = (max-min)*(this->c.at(i-1)-this->c.at(i+1))/(4*i);
    for (; i<=m+1; i++) covi.at(i) = (max-min)*(this->c.at(i-1));
    return ChebyshevApproximation(std::move(covi), min, max);
}

double ChebyshevApproximation::integrate() const {
    double s=0;
    for (int i=1; i<=std::trunc((m+1)/2); i++) s+=(c.at(2*i-2)-c.at(2*i))/(2*i-1);
    return (DajMax()-DajMin())*s/2;
}

double ChebyshevApproximation::derivative(double x) const {
    if (x<min || x>max) throw std::domain_error("Bad argument");
    return this->derivative()(x);
}

double ChebyshevApproximation::integrate(double a, double b) const {
    int t=1;
    if (a>b) {std::swap(a,b);t=-1;}
    if (a<min || b>max) throw std::domain_error("Bad interval");
    ChebyshevApproximation integra(this->antiderivative());
    return (integra(b)-integra(a))*t;
    /*
    double s=0;
    for (int i=1; i<=std::trunc((m+1)/2); i++) s+=(c.at(2*i-2)-c.at(2*i))/(2*i-1);
    return (b-a)*s/2;*/
}


int main() {
    //AT4 - ChebyshevApproximation - 4
const double PI13 = 4 * std::atan(1);
auto funsin3 = [](double x) { return std::sin(x); };
ChebyshevApproximation sinch3(funsin3, 0, PI13, 10);
sinch3.trunc(0.001);
std::cout << funsin3(1) << " " << sinch3(1)<<std::endl;
std::cout <<funsin3(1) << " " << sinch3.derivative(1) << " " << sinch3.derivative()(1)
<< " " << sinch3.derivative().derivative(1)<<std::endl;
std::cout << sinch3.integrate(0, PI13 / 2) << " " << sinch3.integrate();
    return 0;
}