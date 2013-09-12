#ifndef DIFFUSION1DP_H
#define DIFFUSION1DP_H

#include "imc_dfm.h"
#include "tdma.h"

#include "functor1d.h"
#include "boundary.h"


// Classe para coleta de dados dos problemas
class DiffusionData
{
    public:
        tFloat k;
        tFloat ro;
        tFloat C;
        tFloat alpha;
        Functor1D *heatSource;

        DiffusionData(){
            k = 0.0;
            alpha = 0.0;
            ro = 0.0;
            C = 0.0;
            heatSource = new PolynomialConstant(0.0);
        }
};


// Classe para resolução dos problemas de Difusão 1D Transiente
class Diffusion1Dt
{
    public:
        Boundary1D boundaryLeft, boundaryRight;
        DiffusionData data;
        Functor1D *T0;
        tFloat l, hx, ht, tf;
        tInteger nx, nt;
        tFloat **T;

    public:
        TDMA equationsSystem;
        tFloat *Tm, maxValue;

        Diffusion1Dt(){}

        Diffusion1Dt(tInteger nodes, Boundary1D left, Boundary1D right, DiffusionData data, tInteger nt, Functor1D *T0, tFloat tf);

        void plotSolution(tInteger nti, Functor1D &analyticalSolution);
        void printTm(Functor1D &analyticalSolution);
        void plotTm(Functor1D &analyticalSolution);

        void solver(void);

};



#endif // DIFFUSION1DP_H
