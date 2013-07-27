#include "imc_dfm.h"

std::string print(tFloat value)
{
    char str[1000];
    quadmath_snprintf(str, 1000, Q_FORMAT,value);
    return std::string(str);
}

tFloat p_u(tFloat s1, tFloat s2, tFloat s3, tFloat h1, tFloat h2, tFloat h3)
{
    return logq((s2-s3)/(s1-s2))/logq(h3/h2);
}

void GaussSeidel(int n, tFloat *T, tFloat **eq)
{
    T[0] = (eq[0][2]*T[1] + eq[0][2])/eq[0][0];
    for(int i=1; i<n-1; i++)
        T[i] = (eq[i][1]*T[i-1] + eq[i][2]*T[i+1] + eq[i][2])/eq[i][0];
    T[n-1] = (eq[n-1][1]*T[n-2] + eq[n-1][2])/eq[n-1][0];
}

tFloat Residual(int n, tFloat *T, tFloat **eq)
{
    tFloat r = 0.0q;
    r += fabsq(eq[0][2]*T[1] + eq[0][2] - eq[0][0]*T[0]);
    for(int i=1; i<n-1; i++)
        r += fabsq(eq[i][1]*T[i-1] + eq[i][2]*T[i+1] + eq[i][2] - eq[i][0]*T[i]);
    r+= fabsq(eq[n-1][1]*T[n-2] + eq[n-1][2] - eq[n-1][0]*T[n-1]);

    return r;
}

tFloat Ugci(tFloat s1, tFloat s2, tFloat q, tFloat pu, tFloat pl, tFloat Fs)
{
    if(pu<pl)
        return Fs*fabsq(s1-s2)/(powq(q,pu)-1.0q);
    else
        return Fs*fabsq(s1-s2)/(powq(q,pl)-1.0q);
}


