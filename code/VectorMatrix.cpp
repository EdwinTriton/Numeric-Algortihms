#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

/*This file contains implementations of a Vector and Matrix class*/

const char greska[]="Incompatible formats";

class Vector {
    std::vector<double> *elementi;
    public:
    explicit Vector(int n);
    Vector(std::initializer_list<double> l);
    Vector(const Vector& v):elementi(new std::vector<double>(*v.elementi)) {}
    int NElems() const {return elementi->size();}
    double &operator[](int i) {return (*elementi)[i];}
    double operator[](int i) const {return (*elementi)[i];}
    double &operator()(int i) {if (i<=0 || i>NElems()) throw std::range_error("Invalid index");return (*elementi)[i-1];}
    double operator()(int i) const{if (i<=0 || i>NElems()) throw std::range_error("Invalid index");return (*elementi)[i-1];}
    double Norm() const {double x=0;for (int i=0; i<NElems(); i++) x+=(*elementi)[i]*(*elementi)[i];return sqrt(x);}
    friend double VectorNorm(const Vector &v) {return v.Norm();}
    double GetEpsilon() const {return 10*Norm()*std::numeric_limits<double>::epsilon();}
    void Print(char separator = '\n', double eps = -1) const;
    friend void PrintVector(const Vector &v, char separator = '\n',double eps = -1) {v.Print(separator,eps);}
    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v);
    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v);
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s);
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s);
    ~Vector() {delete elementi;}
    Vector& operator=(const Vector& v) {(*elementi)=(*v.elementi);return *this;}
    Vector(Vector&& v):elementi(v.elementi) {v.elementi=nullptr;}
};

Vector::Vector(int n) {
    if (n<=0) throw std::range_error("Bad dimension");
    elementi=new std::vector<double>(n);
}

Vector::Vector(std::initializer_list<double> l) {
    if (l.size()==0) throw std::range_error("Bad dimension");
    elementi=new std::vector<double>(l);
}

void Vector::Print(char separator, double eps) const {
    if (eps<0) eps=GetEpsilon();
    int i;
    for (i=0; i<NElems()-1; i++) {
        double x=(*elementi)[i];
        if (fabs(x)<eps) x=0;
        std::cout<<x<<separator;
    }
    std::cout<<(*elementi)[i];
    if (separator=='\n') std::cout<<'\n';
}

Vector operator +(const Vector &v1, const Vector &v2) {
    if (v1.NElems()!=v2.NElems()) throw std::domain_error(greska);
    Vector v3(v1);
    for (int i=0; i<v1.NElems(); i++) (*v3.elementi)[i]+=(*v2.elementi)[i];
    return v3;
}

Vector operator -(const Vector &v1, const Vector &v2) {
    if (v1.NElems()!=v2.NElems()) throw std::domain_error(greska);
    Vector v3(v1);
    for (int i=0; i<v1.NElems(); i++) (*v3.elementi)[i]-=(*v2.elementi)[i];
    return v3;
}

Vector& Vector::operator +=(const Vector &v) {
    if (NElems()!=v.NElems()) throw std::domain_error(greska);
    for (int i=0; i<NElems(); i++) (*elementi)[i]+=(*v.elementi)[i];
    return *this;
}

Vector& Vector::operator -=(const Vector &v) {
    if (NElems()!=v.NElems()) throw std::domain_error(greska);
    for (int i=0; i<NElems(); i++) (*elementi)[i]-=(*v.elementi)[i];
    return *this;
}

Vector operator *(double s, const Vector &v) {
    Vector izlaz(v);
    for (int i=0; i<izlaz.NElems(); i++) (*izlaz.elementi)[i]*=s;
    return izlaz;
}

Vector operator*(const Vector& v, double s) {
    return s*v;
}

Vector& Vector::operator *=(double s) {
    for (int i=0; i<NElems(); i++) (*elementi)[i]*=s;
    return *this;
}

double operator *(const Vector &v1, const Vector &v2) {
    if (v1.NElems()!=v2.NElems()) throw std::domain_error(greska);
    double x=0;
    for (int i=0; i<v1.NElems(); i++) x+=(*v1.elementi)[i]*(*v2.elementi)[i];
    return x;
}

Vector& Vector::operator /=(double s) {
    if (!s) throw std::domain_error("Division by zero");
    for (int i=0; i<NElems(); i++) (*elementi)[i]/=s;
    return *this;
}

Vector operator /(const Vector &v, double s) {
    Vector k(v);
    k/=s;
    return k;
}

class Matrix {
    std::vector<std::vector<double> > elemente;
    public:
    Matrix(int m, int n) {if (m<=0 || n<=0) throw std::range_error("Bad dimension");elemente=std::vector<std::vector<double> >(m,std::vector<double>(n));}
    Matrix(const Vector &v) {elemente=std::vector<std::vector<double> >(v.NElems(),std::vector<double>(1));for (int i=0; i<v.NElems(); i++) elemente[i][0]=v[i];}
    Matrix(std::initializer_list<std::vector<double>> l);
    int NRows() const {return elemente.size();}
    int NCols() const {if (!elemente.size()) return 0;return elemente.at(0).size();}
    double *operator[](int i) {return elemente[i].data();}
    const double *operator[](int i) const {return elemente[i].data();}
    double &operator()(int i, int j) {if (i<=0 || i>NRows() || j<=0 || j>NCols()) throw std::range_error("Invalid index");return elemente[i-1][j-1];}
    double operator()(int i, int j) const {if (i<=0 || i>NRows() || j<=0 || j>NCols()) throw std::range_error("Invalid index");return elemente[i-1][j-1];}
    double Norm() const ;
    friend double MatrixNorm(const Matrix &m) {return m.Norm();}
    double GetEpsilon() const {return 10*Norm()*std::numeric_limits<double>::epsilon();}
    void Print(int width = 10, double eps = -1) const;
    friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {m.Print(width,eps);}
    /*nesto nesto opasno friend definisat u klasu*/
    friend Matrix operator +(Matrix m1, const Matrix &m2) {return m1+=m2;}
    Matrix &operator +=(const Matrix &m);
    friend Matrix operator -(Matrix m1, const Matrix &m2) {return m1-=m2; }
    Matrix &operator -=(const Matrix &m);
    friend Matrix operator *(double s, const Matrix &m) {return m*s;}
    friend Matrix operator *(Matrix m, double s) {return m*=s;}
    Matrix &operator *=(double s);
    friend Matrix operator *(Matrix m1, const Matrix &m2) {return m1*=m2;}
    Matrix &operator *=(const Matrix &m);
    friend Vector operator *(const Matrix &m, const Vector &v);
    friend Matrix Transpose(Matrix m) {m.Transpose(); return m;}
    void Transpose();
};


double Matrix::Norm() const {
    double x=0;
    for (int i=0; i<NRows(); i++) {
        for (int j=0; j<NCols(); j++) {
            x+=(*this)[i][j]*(*this)[i][j];
        }
    }
    return sqrt(x);
}

void Matrix::Print(int width , double eps) const {
    if (eps<0) eps=GetEpsilon();
    for (int i=0; i<NRows(); i++) {
        for (int j=0; j<NCols(); j++) {
            double x=(*this)[i][j];
            if (fabs(x)<eps) x=0;
            std::cout<<std::right<<std::setw(width)<<x;
        }
        std::cout<<'\n';
    }
}

Matrix& Matrix::operator +=(const Matrix &m) {
    if (NRows()!=m.NRows() || NCols()!=m.NCols()) throw std::domain_error(greska);
    for (int i=0; i<NRows(); i++) {
        for (int j=0; j<NCols(); j++) {
            (*this)[i][j]+=m[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator -=(const Matrix &m) {
    if (NRows()!=m.NRows() || NCols()!=m.NCols()) throw std::domain_error(greska);
    for (int i=0; i<NRows(); i++) {
        for (int j=0; j<NCols(); j++) {
            (*this)[i][j]-=m[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator *=(double s) {
    for (int i=0; i<NRows(); i++) {
        for (int j=0; j<NCols(); j++) {
            (*this)[i][j]*=s;
        }
    }
    return *this;
}

Matrix& Matrix::operator *=(const Matrix &m) {
    if (NCols()!=m.NRows() ) throw std::domain_error(greska);
    Matrix help=(*this);
    for (int i=0; i<NRows(); i++) {
        for (int j=0; j<NCols(); j++) {
            double c=0;
            for (int k=0; k<NCols(); k++) c+=(*this)[i][k]*m[k][j];
            help[i][j]=c;
        }
    }
    elemente=std::move(help.elemente);
    return *this;
}

Vector operator *(const Matrix &m, const Vector &v) {
    if (v.NElems()!=m.NCols()) throw std::domain_error(greska);
    Vector pomocni(m.NRows());
    for (int i=0; i<pomocni.NElems(); i++) for (int j=0; j<m.NCols(); j++) pomocni(i+1)+=m(i+1,j+1)*v(j+1);
    return pomocni;
}

void Matrix::Transpose() {
    if (NRows()==NCols())
    for (int i=0; i<NRows(); i++) for (int j=i+1; j<NCols(); j++) std::swap((*this)[i][j],(*this)[j][i]);
    else {
        std::vector<std::vector<double> > pomocna(NCols(),std::vector<double> (NRows() )  );
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                pomocna[j][i]=(*this)[i][j];
            }
        }
        elemente=std::move(pomocna);
    }
}

Matrix::Matrix(std::initializer_list<std::vector<double>> l) {
    if (!l.size() || !l.begin()->size()) throw std::range_error("Bad dimension");
    int k=(*l.begin()).size(),i=0;
    std::vector<std::vector<double>> pomoc(l.size(),std::vector<double>(k));
    for (const std::vector<double>* p=l.begin(); p!=l.end(); p++,i++) {
        if (p->size()!=k) throw std::logic_error("Bad matrix");
        pomoc[i]=*p;
    }
    elemente=std::move(pomoc);
}

int main()  {
    Vector v1={1,2,3,4,5,6,7,8,9,10},v2={5,4,3,2,1,0,-1,-2,-3,-4},v3(v1);
    Matrix mat1={{1,2,3},{4,5,6},{7,8,9}},mat2={{4,5,6},{9,8,7},{5,5,5}},mat3(mat1),mat4={{1,2,3,4,5},{6,7,8,9,10}},mat5(v1);
    double s=2,norma=v1.Norm(),eps=v1.GetEpsilon();
    /*Testiranje da li se moze grbava matrica formirati*/
    try {
        Matrix {{1,2,3},{1,2}};
        } catch(...) 
        {std::cout<<"Ne moze se formirati grbava matrica\n";}
    try {
        std::cout<<"Prvi el:"<<v1[0]<<' '<<v1(1)<<std::endl;
        /*izuzetak pri koristenje ne dozvoljenog indexa*/
        v1(11);
    }
    catch (...) {
        std::cout<<"Nedopusten index\n";
    }
    if (v1.Norm()!=VectorNorm(v1)) std::cout<<"Ne valja bar jedna od funkcija normi\n";
    if (eps!=10*v1.Norm()*std::numeric_limits<double>::epsilon()) std::cout<<"Greska epsilona v1\n";
    std::cout<<"V1: ";v1.Print(' ');
    std::cout<<"\nV2: ";v2.Print(' ');
    std::cout<<"\nV3: ";v3.Print(' ');
    std::cout<<std::endl;
    (v1+v2).Print(' ');
    std::cout<<std::endl;
    (v3+=v2).Print(' ');
    std::cout<<std::endl;
    (v1-v2).Print(' ');
    std::cout<<std::endl;
    ((v3-=v2)-=v2).Print(' ');
    std::cout<<std::endl;
    (v1*s).Print(' ');
    std::cout<<std::endl;
    ((v3+=v2)*=s).Print(' ');
    std::cout<<std::endl;

    /*Testiranje Matrice:*/
    std::cout<<"Matrica 1:\n";
    mat1.Print();
    std::cout<<"Matrica 2:\n";
    mat2.Print();
    std::cout<<"Matrica 3:\n";
    mat3.Print();
    std::cout<<"Matrica 4:\n";
    mat4.Print();
    std::cout<<"Testiranje matrice: "<<std::endl;
    (mat1+mat2).Print();
    std::cout<<std::endl;
    (mat3+=mat2).Print();
    std::cout<<std::endl;
    (mat1-mat2).Print();
    std::cout<<std::endl;
    ((mat3-=mat2)-=mat2).Print(' ');
    std::cout<<std::endl;
    (mat1*s).Print(' ');
    std::cout<<std::endl;
    ((mat3+=mat2)*=s).Print(' ');
    std::cout<<std::endl;
    try {
        (mat3+mat1).Print();
    }
    catch (...) {
        std::cout<<"Nedozvoljen  operator + za nejednake dimenzije matrice\n";
    }
    try {
        (mat3-mat1).Print();
    }
    catch (...) {
        std::cout<<"Nedozvoljen  operator - za nejednake dimenzije matrice\n";
    }
    try {
        (mat3*mat1).Print();
    }
    catch (...) {
        std::cout<<"Nedozvoljen  operator * za matrice nepravilni dimenzija\n";
    }
    std::cout<<"Matrica 4:\n";
    mat4.Print();
    std::cout<<"Transponovana matrica:\n";
    Transpose(mat4).Print(5);
    std::cout<<std::endl;
    mat4.Transpose();
    mat4.Print(5,2);
    if (MatrixNorm(mat5)!=mat5.Norm() || norma!=mat5.Norm()) std::cout<<"Greska kod norma matrice 5\n";
    if (MatrixNorm(mat3)!=mat3.Norm()) std::cout<<"Greska kod norma matrice 3\n";
    if (mat5.GetEpsilon()!=10*mat5.Norm()*std::numeric_limits<double>::epsilon() || mat5.GetEpsilon()!=eps) std::cout<<"Greska matricnog epsilona matrica 5";
    return 0;
}