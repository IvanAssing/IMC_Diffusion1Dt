#include "functor1d.h"

Functor1D::Functor1D()
{
}

PolynomialConstant::PolynomialConstant(tFloat _a0)
    :a0(_a0)
{

}

tFloat PolynomialConstant::operator ()(tFloat x)
{
    return a0;
}

PolynomialLinear::PolynomialLinear(tFloat _a0, tFloat _a1)
    :a0(_a0), a1(_a1)
{

}

tFloat PolynomialLinear::operator ()(tFloat x)
{
    return a1*x+a0;
}

PolynomialQuadratic::PolynomialQuadratic(tFloat _a0, tFloat _a1, tFloat _a2)
    :a0(_a0), a1(_a1), a2(_a2)
{

}

tFloat PolynomialQuadratic::operator ()(tFloat x)
{
    return (a2*x+a1)*x+a0;
}

PolynomialCubic::PolynomialCubic(tFloat _a0, tFloat _a1, tFloat _a2, tFloat _a3)
    :a0(_a0), a1(_a1), a2(_a2), a3(_a3)
{

}

tFloat PolynomialCubic::operator ()(tFloat x)
{
    return ((a3*x+a2)*x+a1)*x+a0;
}

PolynomialQuartic::PolynomialQuartic(tFloat _a0, tFloat _a1, tFloat _a2, tFloat _a3, tFloat _a4)
    :a0(_a0), a1(_a1), a2(_a2), a3(_a3), a4(_a4)
{

}

tFloat PolynomialQuartic::operator ()(tFloat x)
{
    return (((a4*x+a3)*x+a2)*x+a1)*x+a0;
}


