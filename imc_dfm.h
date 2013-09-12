#ifndef IMC_DFM_H
#define IMC_DFM_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cmath>

extern "C" {
#include "quadmath.h" // Biblioteca para Precisão Quadrúpla
}

#define QUAD_PRECISION // Ativar para habilitar Precisão Quadrúpla

#ifdef QUAD_PRECISION
typedef __float128 tFloat; // Tipo Ponto flutuante com precisão Quadrúpla (128 bits)
#else
typedef double tFloat; // Tipo Ponto flutuante com precisão Dupla (64 bits)
//typedef float tFloat; // Tipo Ponto flutuante com precisão Simples (32 bits)
#endif

typedef long int tInteger; // Tipo Inteiro (64 bits)

#define C_FLOAT(x) static_cast<tFloat>(x)

#ifdef QUAD_PRECISION
#define OUT_FLOAT_PRECISION 33 // Largura do campo de impressão
#else
#define OUT_FLOAT_PRECISION 14
#endif

#define Q_FORMAT "%.33Qe" // Mostrar número em precisão quadrúpla com 33 casas decimais

#define OUT_FLOAT_WIDTH (OUT_FLOAT_PRECISION+10)
#define OUT_TXT 10

#define STR_PAREDE_PLANA "DIFUSÃO DE CALOR EM PAREDE PLANA EM REGIME TRANSIENTE"


tFloat p_u(tFloat s1, tFloat s2, tFloat s3, tFloat h1, tFloat h2, tFloat h3);
void GaussSeidel(int n, tFloat *T, tFloat **eq);
tFloat Residual(int n, tFloat *T, tFloat **eq);
tFloat Ugci(tFloat s1, tFloat s2, tFloat q, tFloat pu, tFloat pl = 2.0q, tFloat Fs = 3.0q);


#define QtoD(value) static_cast<double>(value)

std::string print(tFloat value); // 33 casas
std::string print2(tFloat value); // 3 casas


#endif // IMC_DFM_H
