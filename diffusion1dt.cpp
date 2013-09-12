#include "diffusion1dt.h"

Diffusion1Dt::Diffusion1Dt(tInteger _nodes, Boundary1D _left,
                           Boundary1D _right, DiffusionData _data, tInteger _nt, Functor1D *_T0, tFloat _tf)
{
    boundaryLeft = _left;
    boundaryRight = _right;
    data = _data;
    nx = _nodes;
    nt = _nt+1;
    T0 = _T0;
    tf = _tf;

    // Calcula L e h
    l = boundaryRight.x - boundaryLeft.x;
    hx = l/(nx - 1);
    ht = tf/(nt-1);

    T = new tFloat*[nt];
    for(int i=0; i<nt; i++)
        T[i] = new tFloat[nx];

    Tm = new tFloat[nt];
}


void Diffusion1Dt::solver(void)
{
    // Configura o número de equações do solver TDMA
    equationsSystem.setMaxEquations(nx);

    // DISCRETIZAÇÃO PARA O PROBLEMA DE DIFUSÃO DE CALOR TRANSIENTE EM PAREDA PLANA

    // Aplica a condição de contorno esquerdo
    if(boundaryLeft.type == Dirichlet)
        equationsSystem(1.0q, 0.0q, 0.0q, boundaryLeft.bcValue); // Se Dirichlet

    for(tInteger i=0; i<nx; i++)
        T[0][i] = T0->operator ()(boundaryLeft.x+i*hx);

    // Aplica discretização com CDS-2
    for(tInteger i=1; i<nx-1; i++)
        equationsSystem(1.0q + data.alpha*ht/(hx*hx),
                        data.alpha*ht/(2.0q*hx*hx),
                        data.alpha*ht/(2.0q*hx*hx),
                        0.0q);

    // Aplica a condição de contorno direito
    if(boundaryRight.type == Dirichlet)
        equationsSystem(1.0q, 0.0q, 0.0q, boundaryRight.bcValue); // Se Dirichlet

    // Solver
    for(tInteger ti=1; ti<nt; ti++){
        equationsSystem.T = T[ti];
        // Atualiza o termo fonte
        for(tInteger i=1; i<nx-1; i++)
            equationsSystem.setBp(i, data.alpha*ht/(2.0q*hx*hx)*(T[ti-1][i-1] + T[ti-1][i+1]) +
                    (1.0q - data.alpha*ht/(hx*hx))*T[ti-1][i]);
        // Resolve o sistema de equações lineares
        equationsSystem.solver();
    }

    // Calcula a temperatura média
    for(tInteger ti=0; ti<nt; ti++){
        Tm[ti] = 0.0q;
        for(tInteger i=1; i<nx; i++)
            Tm[ti] += (T[ti][i-1]+T[ti][i]);
        Tm[ti] *= 0.5q*hx/l;
    }

}

void Diffusion1Dt::plotSolution(tInteger nti, Functor1D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "image.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nx; i++)
        file1<<print(boundaryLeft.x+i*hx)<<"\t"<<print(T[nti][i])<<std::endl;
    file1.close();

    // Solução Analítica
    std::ofstream file2(dat2_filename.c_str());
    tFloat hx_10 = l/(10*nx-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*nx; i++){
        tFloat xp = boundaryLeft.x+i*hx_10;
        file2<<static_cast<double>(xp)<<"\t"<<static_cast<double>(analyticalSolution(xp))<<std::endl;
    }
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n";

        file3 << "set title \""<<STR_PAREDE_PLANA<<"\\nResolução com MDF \"\n"
                 "set xlabel 'x'\n"
                 "set ylabel 'T(x, "<<QtoD(ht*nti)<<")'\n";

    file3 <<

             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}


// Imprime os resultados
void Diffusion1Dt::printTm(Functor1D &analyticalSolution)
{
        for(tInteger ti=0; ti<nt; ti++)
            std::cout<<std::endl<<ti<<"\t"<<print2(ti*ht)<<"\t"<<print(Tm[ti])<<"\t"<<
                       print(analyticalSolution(ti*ht))<<"\t"<<print(analyticalSolution(ti*ht)-Tm[ti]);

}

void Diffusion1Dt::plotTm(Functor1D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "image.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nt; i++)
        file1<<print(i*ht)<<"\t"<<print(Tm[i])<<std::endl;
    file1.close();

    // Solução Analítica
    std::ofstream file2(dat2_filename.c_str());
    tFloat h_10 = tf/(10*nt-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*nt; i++){
        file2<<static_cast<double>(i*h_10)<<"\t"<<static_cast<double>(analyticalSolution(i*h_10))<<std::endl;
    }
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n";

        file3 << "set title \""<<STR_PAREDE_PLANA<<"\\nResolução com MDF \"\n"
                 "set xlabel 't'\n"
                 "set ylabel 'Temperatura Média(t)'\n";

    file3 <<

             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}

