#ifndef FUNCTOR1D_H
#define FUNCTOR1D_H

#include "imc_dfm.h"

// Objeto Funcão (Super Classe Abstrata)
class Functor1D
{
    public:
        Functor1D();
        virtual tFloat operator()(tFloat x) = 0;
};

// Objeto Polinômio Constante
class PolynomialConstant : public Functor1D
{
    public:
        tFloat a0;
        PolynomialConstant(tFloat a0);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 1
class PolynomialLinear : public Functor1D
{
    public:
        tFloat a0, a1;
        PolynomialLinear(tFloat a0, tFloat a1);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 2
class PolynomialQuadratic : public Functor1D
{
    public:
        tFloat a0, a1, a2;
        PolynomialQuadratic(tFloat a0, tFloat a1, tFloat a2);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 3
class PolynomialCubic : public Functor1D
{
    public:
        tFloat a0, a1, a2, a3;
        PolynomialCubic(tFloat a0, tFloat a1, tFloat a2, tFloat a3);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 4
class PolynomialQuartic : public Functor1D
{
    public:
        tFloat a0, a1, a2, a3, a4;
        PolynomialQuartic(tFloat a0, tFloat a1, tFloat a2, tFloat a3, tFloat a4);
        virtual tFloat operator()(tFloat x);
};



class Sine : public Functor1D
{
    public:
        tFloat a, b, c;

        Sine(tFloat _a, tFloat _b, tFloat _c):a(_a), b(_b), c(_c){}

        tFloat operator()(tFloat x)
        {
            return a*sinq(b*x+c);
        }
};

class SpecialFunctorTxt : public Functor1D
{
    public:
        tFloat a, b, c;

        SpecialFunctorTxt(tFloat _a, tFloat _b, tFloat _c):a(_a), b(_b), c(_c){}

        tFloat operator()(tFloat x)
        {
            return a*sinq(b*x)*expq(c);
        }
};

class SpecialFunctorTm : public Functor1D
{
    public:
        tFloat a, b;

        SpecialFunctorTm(tFloat _a, tFloat _b):a(_a), b(_b){}

        tFloat operator()(tFloat x)
        {
            return a*expq(b*x);
        }
};

#endif // FUNCTOR1D_H
