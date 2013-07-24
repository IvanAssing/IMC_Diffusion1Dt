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

    /*
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
*/

    /*
    // >[QUESTÃO 4.3]
    int n[3] = {11, 21, 41};
    //int n[3] = {1000001, 2000001, 4000001};
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

    Diffusion1Dp *mesh = new Diffusion1Dp[3];

    for(int i=0; i<3; i++){
        mesh[i] = Diffusion1Dp(ParedePlana, n[i], left, right, data);
        mesh[i].solver();
    }


    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
    std::cout<<std::endl<<std::setw(5)<<std::right<<"N";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"T(1/2)";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"Tm";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"q\"(0)";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"q\"(L)";
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');


    for(int i=0; i<3; i++){
        std::cout<<std::endl<<std::setw(5)<<n[i];
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[i].equationsSystem.getT(n[i]/2 + 1));
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[i].averageValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[i].heatFlowLeft);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[i].heatFlowRight);

    }

    std::cout<<std::endl<<std::setw(5)<<"q";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[0].get_h()/mesh[1].get_h());
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[0].get_h()/mesh[1].get_h());
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[0].get_h()/mesh[1].get_h());
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(mesh[0].get_h()/mesh[1].get_h());

    std::cout<<std::endl<<std::setw(5)<<"pu";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(p_u(mesh[2].equationsSystem.getT(n[2]/2 + 1),
            mesh[1].equationsSystem.getT(n[1]/2 + 1), mesh[0].equationsSystem.getT(n[0]/2 + 1),
            mesh[2].get_h(), mesh[1].get_h(), mesh[0].get_h()));
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(p_u(mesh[2].averageValue, mesh[1].averageValue, mesh[0].averageValue,
            mesh[2].get_h(), mesh[1].get_h(), mesh[0].get_h()));
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(p_u(mesh[2].heatFlowLeft, mesh[1].heatFlowLeft, mesh[0].heatFlowLeft,
            mesh[2].get_h(), mesh[1].get_h(), mesh[0].get_h()));
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(p_u(mesh[2].heatFlowRight, mesh[1].heatFlowRight, mesh[0].heatFlowRight,
            mesh[2].get_h(), mesh[1].get_h(), mesh[0].get_h()));


    std::cout<<std::endl<<std::setw(5)<<"erro";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(AS_Temperature_PP(0.5q) - mesh[2].equationsSystem.getT(n[2]/2 + 1));
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(AS_tm - mesh[2].averageValue);
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(AS_q0 - mesh[2].heatFlowLeft);
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<print(AS_ql - mesh[2].heatFlowRight);



    //mesh.printSolution(AS_Temperature_PP);
    //mesh.printSecondaryResults(AS_tm, AS_q0, AS_ql, 0.0q, 1.0q);

    //mesh.plotSolution(AS_Temperature_PP);

    // <[QUESTÃO 4.3]

    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
    std::cout<<"\n\n\n";
    */

    /*
    // >[QUESTÃO 4.4]
    tFloat s0 = -0.5q;
    tFloat s1 = -1.5q;
    tFloat s2 = -1.0q;
    DiffusionData data;
    data.k = 1.0q;
    data.heatSource = new PolynomialQuadratic(s0, s1, s2); // s0 + s1*x + s2*x2

    Boundary1D left(0.0q, Dirichlet, 0.0q); // x=0, Condição de Contorno de Dirichlet, T(0) = 0
    Boundary1D right(1.0q, Dirichlet, 1.0q); // x=1, Condição de Contorno de Dirichlet, T(1) = 1

    int n = 21;

    Diffusion1Dp mesh(ParedePlana, n, left, right, data);

    mesh.solver();

    tFloat *T = new tFloat[n];

    // Estimativa para solução (reta)
    tFloat _m = (right.bcValue-left.bcValue)/(right.x-left.x);
    tFloat _h = (right.x-left.x)/(n-1);

    for(int i=0; i<n; i++)
        T[i] = left.bcValue + _m*(i*_h+left.x);

    int itmax = 1000000, it=0;
    tFloat *L = new tFloat[itmax];
    tFloat itol = 1.0e-35q;

    L[0] = Residual(n, T, mesh.equationsSystem.equation);

    //std::cout<<"\n"<<0<<"\t"<<static_cast<double>(L[0])<<"\t"<<static_cast<double>(itol);

    do{
        ++it;
        GaussSeidel(n, T, mesh.equationsSystem.equation);
        L[it] = Residual(n, T, mesh.equationsSystem.equation);

        //std::cout<<"\n"<<it<<"\t"<<static_cast<double>(L[it]);

    }while(it<itmax && L[it]>itol);


    const std::string cmd_filename = "plotconfig_it.gnu";
    const std::string pic_filename = "iteractive_it.png";
    const std::string dat1_filename = "data_it.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=1; i<it; i++)
        file1<<i<<"\t"<<static_cast<double>(L[i]/L[0])<<std::endl;
    file1.close();


    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set logscale y\n"
             "set format y \"10^{%L}\" \n"
             "set lmargin 10 \n";
    file3 << "set title \""<<STR_PAREDE_PLANA<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set ylabel \"L^{n}/L^{0}\" \n"
             "set xlabel 'Número de iterações'\n";

    file3 <<"plot '" <<dat1_filename<<"' t\"\" with lines lt 2 lc 1 lw 1";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());


    // <[QUESTÃO 4.4]
*/

    // >[QUESTÃO 4.5]
    DiffusionData data;
    data.k = 1.0q;
    data.heatSource = new PolynomialConstant(0.0q); // s0 + s1*x + s2*x2

    PolynomialLinear AS_T(0.0q, 1.0q); // solução analítica

    Boundary1D left(0.0q, Dirichlet, 0.0q); // x=0, Condição de Contorno de Dirichlet, T(0) = 0
    Boundary1D right(1.0q, Dirichlet, 1.0q); // x=1, Condição de Contorno de Dirichlet, T(1) = 1

    int nmax = 10000000;

    int itmax = log10(nmax);

    std::cout<<"\n"<<itmax<<"\t"<<nmax;

    tFloat *erro = new tFloat[itmax];

    int n =10;

    for(int i=0; i<itmax; i++, n*=10){

        Diffusion1Dp mesh(ParedePlana, n+1, left, right, data);

        mesh.solver();

        erro[i] = fabsq(0.5q - mesh.equationsSystem.getT(n/2+1));

        std::cout<<"\n"<<i<<"\t"<<static_cast<double>(erro[i])<<"\t"<<n+1;

    }


    const std::string cmd_filename = "plotconfig_45.gnu";
    const std::string pic_filename = "iteractive_45.png";
    const std::string dat1_filename = "data_45.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=1; i<itmax; i++)
        file1<<1.0/pow(10., i)<<"\t"<<static_cast<double>(erro[i])<<std::endl;
    file1.close();


    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set logscale xy\n"
             "set format y \"10^{%L}\" \n"
             "set format x \"10^{%L}\" \n"
             "set lmargin 10 \n";
    file3 << "set title \""<<STR_PAREDE_PLANA<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set ylabel \"L^{n}/L^{0}\" \n"
             "set xlabel 'Número de iterações'\n";

    file3 <<"plot '" <<dat1_filename<<"' t\"\" with linespoints lt 2 lc 4 lw 2";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());


    // <[QUESTÃO 4.5]
}

