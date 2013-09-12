#include <iostream>

#include "tdma.h"
#include "functor1d.h"
#include "diffusion1dt.h"

int main()
{
    // >[QUESTÃO 9.1] - DIFUSÃO DE CALOR EM PAREDE PLANA EM REGIME TRANSIENTE
    DiffusionData data;
    data.alpha = 1.0q;

    tFloat T_0 = 1.0q;
    tFloat L = 1.0q;
    tFloat tf = 0.1q;

    Sine T0(T_0, M_PIq/L, 0.0q); // Temperatura inicial

    Boundary1D left(0.0q, Dirichlet, 0.0q); // x=0, Condição de Contorno de Dirichlet, T(0,t) = 0
    Boundary1D right(L, Dirichlet, 0.0q); // x=1, Condição de Contorno de Dirichlet, T(1,t) = 0

    Diffusion1Dt mesh(17, left, right, data, 40, &T0, tf); // Malha

    mesh.solver();

    // Soluções analíticas
    SpecialFunctorTxt AS_T(T_0, M_PIq/L, -data.alpha*M_PIq*M_PIq*tf/(L*L));
    SpecialFunctorTm AS_Tm(2.0q/M_PIq*T_0, -data.alpha*M_PIq*M_PIq/(L*L));

    mesh.plotSolution(40, AS_T);

    mesh.printTm(AS_Tm);
    mesh.plotTm(AS_Tm);

    // Calculo da estimativa do erro de discretização
    Diffusion1Dt meshx1(40, left, right, data, 40, &T0, tf);
    Diffusion1Dt meshx2(20, left, right, data, 40, &T0, tf);
    Diffusion1Dt meshx3(10, left, right, data, 40, &T0, tf);

    meshx1.solver();meshx2.solver();meshx3.solver();
    tFloat Tx1 = meshx1.Tm[40];
    tFloat Tx2 = meshx2.Tm[40];
    tFloat Tx3 = meshx3.Tm[40];

    tFloat pux = p_u(Tx1, Tx2, Tx3, meshx1.hx, meshx2.hx, meshx3.hx);

    tFloat ugcix = Ugci(Tx1, Tx2, meshx2.hx/meshx1.hx, pux);

    Diffusion1Dt mesht1(17, left, right, data, 40, &T0, tf);
    Diffusion1Dt mesht2(17, left, right, data, 20, &T0, tf);
    Diffusion1Dt mesht3(17, left, right, data, 10, &T0, tf);

    mesht1.solver();mesht2.solver();mesht3.solver();
    tFloat Tt1 = mesht1.Tm[40];
    tFloat Tt2 = mesht2.Tm[20];
    tFloat Tt3 = mesht3.Tm[10];

    tFloat put = p_u(Tt1, Tt2, Tt3, mesht1.ht, mesht2.ht, mesht3.ht);

    tFloat ugcit = Ugci(Tt1, Tt2, mesht2.ht/mesht1.ht, put);

    std::cout<<"\n\nTm("<<QtoD(tf)<<") = "<<print(ugcit>ugcix?Tt1:Tx1)<<" +- "<<print(ugcit>ugcix?ugcit:ugcix);

}

