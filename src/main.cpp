#include "../include/C_prb_RR_HMS.h"    

#include<fstream>
#include<iostream>
#include<iomanip>
#include<chrono>
using namespace std;

int main(){
    // Hybrid Multi Shooting
    test_RR_HMS();

    // VectorXd YY(6), YY_ODE(6), xf_ODE(6); 
    // VectorXd x0rv(6), xf(6), x0(42); x0rv << 3., 0., 0., 0., 2.5, 1.;
    // MatrixXd identityMatrix = MatrixXd::Identity(6, 6);
    // MatrixXd PHIf(6, 6), PHIf_ODE(6, 6);
    // x0 << x0rv, identityMatrix.reshaped();
    // auto start_time = chrono::high_resolution_clock::now();
    // C_simulazione sim(42, 1e3, 1, 1e-14, 1e-14);
    // sim.SetY0_Sim(x0);
    // EoM_Kepler_STM EoM(1., 2*M_PI);
    // sim.Start_Sim(EoM);
    // YY = sim.Y.rightCols(1);
    // xf = YY.head(6);
    // PHIf = YY.tail(36).reshaped(6, 6);
    // cout << "xf_RK4 = " << xf.transpose() << endl;
    // cout << "PHIf_RK4 = " << endl << PHIf.transpose() << endl;
    // auto end_time = chrono::high_resolution_clock::now();
    // auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    // double milliseconds_TK4 = static_cast<double>(duration.count())/1e3;

    // start_time = std::chrono::high_resolution_clock::now();
    // C_simulazione_ODE sim_ODE(42, 2, 1, 1e-14, 1e-14);
    // sim_ODE.SetY0_Sim(x0);
    // EoM_Kepler_ODE_STM EoM_ODE(1., 2*M_PI);
    // sim_ODE.Start_Sim(EoM_ODE);
    // YY_ODE = sim_ODE.Y.rightCols(1);
    // xf_ODE = YY_ODE.head(6);
    // PHIf_ODE = YY_ODE.tail(36).reshaped(6, 6);
    // cout << "YY_RK8 = " << xf_ODE.transpose() << endl;
    // cout << "PHIf_RK8 = " << endl << PHIf_ODE << endl;
    // end_time = chrono::high_resolution_clock::now();
    // duration = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // double milliseconds_TK8 = static_cast<double>(duration.count())/1e3;
    // cout << fixed << setprecision(5) << "T_RK4: " << milliseconds_TK4 << " ms" << std::endl;
    // cout << fixed << setprecision(5) << "T_RK8: " << milliseconds_TK8 << " ms" << std::endl;
   
    // auto start_time = chrono::high_resolution_clock::now();
    // Vector3d r0, v0;
    // r0 = x0rv.head(3); v0 = x0rv.tail(3);
    // MatrixXd Qd_k(6, 6), P0(6, 6), Pf(6, 6);
    // Qd_k.setZero(); MatrixXd::Zero(6, 6);

    // double sigmar = 1e-4;
    // double sigmav = 1e-3;

    // P0.setZero();
    // P0.block<3, 3>(0, 0) = sigmar * sigmar * Matrix3d::Identity();
    // P0.block<3, 3>(3, 3) = sigmav * sigmav * Matrix3d::Identity();
    // auto end_time = chrono::high_resolution_clock::now();
    // auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    // double milliseconds = static_cast<double>(duration.count())/1e3;

    // Pf_STM(r0, v0, P0, 2*M_PI, 1, MatrixXd::Zero(6, 6), Pf);
    // cout << "Pf = " << endl << Pf << endl;
    // cout << fixed << setprecision(5) << "Time: " << milliseconds << " ms" << std::endl;
return 0;

}



