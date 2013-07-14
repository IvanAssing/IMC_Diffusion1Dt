#include <iostream>

#include "tdma.h"
#include "functor1d.h"
#include "diffusion1dp.h"

int main()
{
/*
    // >[QUESTÃO 3.1] - SISTEMA DE EQUAÇÕES COM MÉTODO TDMA
    tInteger nEquations = 11; // Número de Equações
    TDMA equationsSystem(nEquations);

    equationsSystem.addEquation(1.0q, 0.0q, 0.0q, 0.0q);

    for(tInteger i=1; i<nEquations-1; i++)
        equationsSystem.addEquation(2.0q, 1.0q, 1.0q, 0.0q); // ap, aw, ae, bp

    equationsSystem.addEquation(1.0q, 0.0q, 0.0q, 1.0q);

    equationsSystem.solver();

    equationsSystem.printCoefficients();
    equationsSystem.printSolution();
    // <[QUESTÃO 3.1]
*/

/*
    // >[QUESTÃO 3.2] - DIFUSÃO DE CALOR EM PAREDE PLANA COM TEMPERATURA PRESCRITA
    tInteger n = 1000001;
    tFloat s0 = -0.5q;
    tFloat s1 = -1.5q;
    tFloat s2 = -1.0q;
    DiffusionData data;
    data.k = 1.0q;
    data.heatSource = new PolynomialQuadratic(s0, s1, s2); // s0 + s1*x + s2*x2

    PolynomialQuartic AS_Temperature_PP(0.0q, 1.0q - s0/2.0q - s1/6.0q - s2/12.0q,
                                        s0/2.0q, s1/6.0q, s2/12.0q); // solução analítica
    tFloat AS_tm = 0.5q - s0/12.0q - s1/24.0q - s2/40.0q; // temperatura média
    tFloat AS_q0 = -1.0q + s0/2.0q + s1/6.0q + s2/12.0q; // fluxo de calor em x=0
    tFloat AS_ql = -(1.0q + s0/2.0q + s1/3.0q + s2/4.0q); // fluxo de calor em x=1

    Boundary1D left(0.0q, Dirichlet, 0.0q); // x=0, Condição de Contorno de Dirichlet, T(0) = 0
    Boundary1D right(1.0q, Dirichlet, 1.0q); // x=1, Condição de Contorno de Dirichlet, T(1) = 1

    Diffusion1Dp mesh(ParedePlana, n, left, right, data);

    mesh.solver();

    //mesh.printSolution(AS_Temperature_PP);
    mesh.printSecondaryResults(AS_tm, AS_q0, AS_ql, 0.0q, 1.0q);

    //mesh.plotSolution(AS_Temperature_PP);
    // <[QUESTÃO 3.2]
*/

/*
    // >[QUESTÃO 3.3] - DIFUSÃO DE CALOR EM PAREDE PLANA COM FLUXO DE CALOR PRESCRITO
    tInteger n = 21;
    tFloat s0 = -0.5q;
    tFloat s1 = -3.0q;
    tFloat s2 = -1.0q;
    DiffusionData data;
    data.k = 1.0q;
    data.heatSource = new PolynomialQuadratic(s0, s1, s2);
    tFloat ql = 0.5q;

    PolynomialQuartic AS_Temperature_PP(0.0q, -(ql/data.k + s0 + s1/2.0q + s2/3.0q),
                                        s0/2.0q, s1/6.0q, s2/12.0q);
    tFloat AS_tm = -0.5*(ql/data.k + s0 + s1/2.0q + s2/3.0q) + s0/6.0q + s1/24.0q + s2/60.0q;
    tFloat AS_q0 = data.k*(ql/data.k + s0 + s1/2.0q + s2/3.0q);
    tFloat AS_ql = ql;

    Boundary1D left(0.0q, Dirichlet, 0.0q); // x=0, Condição de Contorno de Dirichlet, T(0) = 0
    Boundary1D right(1.0q, Neumann, ql); // x=1, Condição de Contorno de Neumann, q"(1) = 1/2

    Diffusion1Dp mesh(ParedePlana, n, left, right, data);

    mesh.solver();

    //mesh.printSolution(AS_Temperature_PP);
    mesh.printSecondaryResults(AS_tm, AS_q0, AS_ql, 0.0q, AS_Temperature_PP(1.0q));

    mesh.plotSolution(AS_Temperature_PP);

    // <[QUESTÃO 3.3]
*/

/*
    // >[QUESTÃO 3.4] - DIFUSÃO DE CALOR EM ALETA
        tInteger n = 1000001;
        DiffusionData data;

        tFloat D = 0.005q;
        data.P = M_PIq*D;
        data.Tb = 100.0q;
        data.H = 100.0q;
        data.Ab = M_PIq*D*D/4.0q;
        data.Tinf = 25.0q;
        data.k = 398.0q;
        tFloat L = 0.2q;

        Boundary1D left(0.0q, Dirichlet, 0.0q);
        Boundary1D right(L, Robin);

        SpecialFunctorTA AS_Temperature_Aleta(data.k, data.H, data.Tinf, data.Tb, data.Ab, data.P, L);
        SpecialFunctorQA AS_HeatFlow_Aleta(data.k, data.H, data.Tinf, data.Tb, data.Ab, data.P, L);

        Diffusion1Dp mesh(Aleta, n, left, right, data);

        mesh.solver();

        //mesh.printSolution(AS_Temperature_Aleta);
        mesh.printSecondaryResults(C_FLOAT(151./240.), AS_HeatFlow_Aleta(0.0q), AS_HeatFlow_Aleta(0.2q));

        //mesh.plotSolution(AS_Temperature_Aleta);

    // <[QUESTÃO 3.4]
*/


    // >[QUESTÃO 3.5] - DIFUSÃO DE QUANTIDADE DE MOVIMENTO LINEAR
    tInteger n = 11;
    DiffusionData data;

    tFloat D = 5.e-2q;
    data.mi = 1.e-3q;
    data.C = -9.6q;

    PolynomialQuadratic AS_u_qml(0.0q, -0.5q*data.C*D/data.mi, 0.5q*data.C/data.mi);
    tFloat AS_u_m = -data.C*D*D/(data.mi*12.0q);
    tFloat AS_u_max = 1.5*AS_u_m;

    Boundary1D left(0.0q); // x=0, CC_default = Dirichlet
    Boundary1D right(D); // x=1, CC_default = Dirichlet

    Diffusion1Dp mesh(QML, n, left, right, data);


    mesh.solver();

    mesh.printSolution(AS_u_qml);
    mesh.printSecondaryResults(AS_u_m, AS_u_max);


    mesh.plotSolution(AS_u_qml);
    // <[QUESTÃO 3.5]


}

