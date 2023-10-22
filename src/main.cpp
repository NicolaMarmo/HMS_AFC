#include "../include/C_prb_RR_HMS.h"    

#include<fstream>
#include<iostream>
#include<iomanip>
#include <chrono>
using namespace std;

int main()
{
    // Hybrid Multi Shooting
    test_RR_HMS();

    // VectorXd YY(6), YY_ODE(6);
    // auto start_time = std::chrono::high_resolution_clock::now();
    // VectorXd x0(6); x0 << 1., 0., 0., 0., 1., 0.;
    // C_simulazione sim(6, 1e3, 1, 1e-14, 1e-14);
    // sim.SetY0_Sim(x0);
    // EoM_Kepler_STM EoM(1., 2*M_PI);
    // sim.Start_Sim(EoM);
    // YY = sim.Y.rightCols(1);
    // cout << "YY_RK4 = " << YY.transpose() << endl;
    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // double milliseconds_TK4 = static_cast<double>(duration.count()) / 1000.0;

    // start_time = std::chrono::high_resolution_clock::now();
    // C_simulazione_ODE sim_ODE(6, 2, 1, 1e-14, 1e-14);
    // sim_ODE.SetY0_Sim(x0);
    // EoM_Kepler_ODE EoM_ODE(1., 2*M_PI);
    // sim_ODE.Start_Sim(EoM_ODE);
    // YY_ODE = sim_ODE.Y.rightCols(1);
    // cout << "YY_RK8 = " << YY_ODE.transpose() << endl;
    // end_time = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // double milliseconds_TK8 = static_cast<double>(duration.count()) / 1000.0;
    // std::cout << std::fixed << std::setprecision(2) << "T_RK4: " << milliseconds_TK4 << " ms" << std::endl;
    // std::cout << std::fixed << std::setprecision(3) << "T_RK8: " << milliseconds_TK8 << " ms" << std::endl;

    return 0;
}



