#include "../include/ut.h"

void propagate_UT(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k, int Nsub,
                        Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1){
    // Q_k è il disturbo discreto addizionale
    
    int n = 6;    //state dimension in 3d

    VectorXd s0(6); s0 << r0_mean, v0_mean;

    C_UT ut(6);

    // genero sigma_points XX(xMean, Pxx)
    MatrixXd XX;
    ut.eval_sigmaPoints(s0, P0, XX);

    // propago: YY = propagate_kep(XX,dt)   
    MatrixXd YY(n, 2*n+1);
    C_simulazione sim(n, 100, 1, 1e-6, 1e-8);
    EoM_CR3BP EoM(mu, dt_sec);
    for (int j=0; j<1+2*n; j++){
        s0 = XX.col(j);

        sim.SetY0_Sim(s0);
        sim.Start_Sim(EoM);
        YY.col(j) = sim.Y.rightCols(1);
    } 


    // ricampiono indietro    
    VectorXd s1_mean;
    ut.inverse_sigmaPoints(YY, s1_mean, P1);

    r1_mean = s1_mean.segment(0,3);
    v1_mean = s1_mean.segment(3,3);

    // aggiungi rumore di processo (Qd)
    P1 = P1 + Qd_k;

    // DEBUG
    // cout << "--------------------------------------" << endl;
    // cout << "XX = \n" << XX << endl;
    // cout << "P1_plus = \n" << P1_plus << endl;
}
