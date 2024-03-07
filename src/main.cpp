#include "../include/C_prb_RR_HMS.h"    

#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace std;

int main(){
    // Hybrid Multi Shooting
    test_RR_HMS();
    // MatrixXd P0(6, 6), Pf(6, 6);
    // Vector3d v_infm, v_infp, r_infm, r_infp;
    // double ToF, DT;
    // r_infm << 640066.5250193276, 647387.6090269776, 185066.4372003882;
    // v_infm << -1.6620084020, -1.6289505750, -0.6192816393;
    // v_infp << -1.84769317509121,    1.45832611672181,	 0.841989674781640;
    // P0.diagonal() << 2.5e7/1.5095e+08, 2.5e7/1.5095e+08, 2.5e7/1.5095e+08, 1e-6/29.6510, 1e-6/29.6510, 1e-6/29.6510;
    // // Propagate_P_FB(P0, v_infm, v_infp, 929000, 3.986004418e5, Pf);
    // ToF_Hyperbola(r_infm, v_infm, 3.986004418e5, ToF);
    // propagateKEP_U(r_infm/1.5095e+08, v_infm/29.6510, 6.637443485871841e+05/5.0909e+06, 3.986004418e5/132712440018, r_infp, v_infp);
    // Pf_STM(r_infm/1.5095e+08, v_infm/29.6510, P0, 0, 3.986004418e5/132712440018, MatrixXd::Zero(6, 6), Pf);
    // cout << "P0 = " << endl << P0 << endl;
    // cout << "Pf = " << endl << Pf << endl;

    // Vector3d r0, v0;
    // MatrixXd P0(6, 6), Pf(6, 6);
    // double dt_sec = 6.637443485871841e+05;
    // r0 << -778347.923501835,	-439171.761750917,	-253563.372835403;
    // v0 << 2.16500000000000,	1.08250000000000,	0.625000000000000;
    // P0.diagonal() << 2.5e7, 2.5e7, 2.5e7, 1e-6, 1e-6, 1e-6;
    // Pf_STM(r0, v0, P0, dt_sec, 3.986004418e5, MatrixXd::Zero(6, 6),
    //         Pf);
    // cout << "Pf = " << endl << Pf << endl;

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



