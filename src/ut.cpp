#include "../include/ut.h"


// void ut_propagation(int n)
// {

//     MatrixXd P_k(n,n);


//     VectorXd x0(3); x0 << 0., 0., 0.;

//     MatrixXd XX_k(n, 2*n+1);
//     VectorXd WW(2*n+1);


//     MatrixXd P0(3,3);

//     // P = L*D*L^T

//     Eigen::MatrixXd A(3,3), P(3,3);
//     // P << 6, 0, 0, 0, 4, 0, 0, 0, 7;
//     A << 1,2,3, 
//          6,3,4, 
//          3,4,5;
//     P = A.transpose()*A;
//     Eigen::MatrixXd L( P.llt().matrixL() );
//     // std::cout << L.col(0) << std::endl;
    
//     // cout << "D = \n" >> P.llt().matrixD() << endl;
    
//     // cout << "P = \n" << P << endl;
//     // cout << "L = \n" << L << endl;
//     // cout << "P-L*L^T = \n" << P-L*L.transpose() << endl;

//     MatrixXd XX(n,1+2*n);

//     XX.col(0)=x0;
//     for (int j=0; j<n; j++)
//         XX.col(j) = x0 + L.col(j);

//     for (int j=0; j<n; j++)
//         XX.col(j+n) = x0 - L.col(j);
    
//     cout << XX << endl;
// }


void propagate_kepler_UT(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k,
                        Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1)
{
    // Q_k è il disturbo discreto addizionale
    
    int n=6;    //state dimension in 3d

    VectorXd s0(6); s0 << r0_mean, v0_mean;

    C_UT ut(6);

    // genero sigma_points XX(xMean, Pxx)
    MatrixXd XX;
    ut.eval_sigmaPoints(s0, P0, XX);


    // propago: YY = propagate_kep(XX,dt)   
    MatrixXd YY(n,2*n+1);
    for (int j=0; j<1+2*n; j++)
    {
        s0 = XX.col(j);
        Vector3d r0 = s0.segment(0,3);
        Vector3d v0 = s0.segment(3,3);
        Vector3d r1, v1;
        propagateKEP_U(r0, v0, dt_sec, mu, r1, v1);
        VectorXd s1(6); s1 << r1,v1; 
        YY.col(j) = s1;
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

void propagate_kepler_UT_MS(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k,
                        Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1)
{
    // Q_k è il disturbo discreto addizionale
    int npar = 50;

    int n=6;    //state dimension in 3d

    VectorXd s0(6); s0 << r0_mean, v0_mean;

    C_UT ut(6);

    VectorXd s1_mean(6);
    VectorXd s1(6);
    Vector3d r0, v0, r1, v1;

    MatrixXd YY(n,2*n+1);
    
    // genero sigma_points XX(xMean, Pxx)
    MatrixXd XX;

    for (int k=0; k<npar; k++)
    {
        ut.eval_sigmaPoints(s0, P0, XX);

        // propago: YY = propagate_kep(XX,dt)   
        for (int j=0; j<1+2*n; j++)
        {
            s0 = XX.col(j);
            r0 = s0.segment(0,3);
            v0 = s0.segment(3,3);
            propagateKEP_U(r0, v0, dt_sec/npar, mu, r1, v1);
            s1 << r1,v1; 
            YY.col(j) = s1;
        } 

        // ricampiono indietro    
        ut.inverse_sigmaPoints(YY, s1_mean, P1);

        s0 = s1_mean;
        P0 = P1;
    }
    r1_mean = s0.segment(0,3);
    v1_mean = s0.segment(3,3);

    // aggiungi rumore di processo (Qd)
    P1 = P0 + Qd_k;

    // DEBUG
    // cout << "--------------------------------------" << endl;
    // cout << "XX = \n" << XX << endl;
    // cout << "P1_plus = \n" << P1_plus << endl;
}

inline MatrixXd GGravity(const Vector3d& r, const double amu)
{
    double x=r(0);
    double y=r(1);
    double z=r(2);

    double r_norm = r.norm();
    double ur5 = 1/pow(r_norm,5);
    double r_norm2 = pow(r_norm,2);
    MatrixXd G(3,3);
    G << 3*pow(x,2)-r_norm2, 3*x*y, 3*x*z,
         3*x*y,   3*pow(y,2)-r_norm2, 3*y*z,
         3*x*z, 3*y*z, 3*pow(z,2)-r_norm2;
        
    G *= amu*ur5;
    return G;
}

//TODO: cercare reference della propagazione STM
/***
 * 
 *      Qd_k = dt *  G*G^T
 * 
 *      nSub =numero di divisioni interne nella propagazione
 **/
void propagate_kepler_STM(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double amu, MatrixXd Qd_k,
    Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1, int nSub)
{

    const int n=6;

    // MatrixXd STM;
    // propagateKEP_U(r0_mean, v0_mean, dt_sec, mu, r1_mean, v1_mean, STM);
    // P1 = STM * P0 * STM.transpose() + Qd_k;

    // if (nSub<1)
    //     nSub=1;


    // //------------------------------------------//
    // propagateKEP_U(r0_mean, v0_mean, dt_sec, amu, r1_mean, v1_mean);
   
    // MatrixXd G0 = GGravity(r0_mean, amu);
    // MatrixXd G = GGravity(r1_mean, amu);
    // MatrixXd PHI_rr = MatrixXd::Identity(3,3) + (2*G0+G)*pow(dt_sec,2)/6.;
    // MatrixXd PHI_rv = MatrixXd::Identity(3,3)*dt_sec +(G0+G)*pow(dt_sec,3)/12.;
    // MatrixXd PHI_vr = (G0+G)*dt_sec/2.;
    // MatrixXd PHI_vv = MatrixXd::Identity(3,3) + (2*G0+G)*pow(dt_sec,2)/6.;

    // MatrixXd STM(6,6); STM << PHI_rr, PHI_rv, PHI_vr, PHI_vv;
    // P1 = STM * P0 * STM.transpose() + Qd_k;
    // //------------------------------------------//

    /* sub-segments */
    Vector3d _r0 = r0_mean;
    Vector3d _v0 = v0_mean;
    MatrixXd _P0 =P0;
    double dt_seg = dt_sec/double(nSub);
    for (int k=0; k<nSub; k++)
    {

        propagateKEP_U(_r0, _v0, dt_seg, amu, r1_mean, v1_mean);
   
        MatrixXd G0 = GGravity(r0_mean, amu);
        MatrixXd G = GGravity(r1_mean, amu);
        MatrixXd PHI_rr = MatrixXd::Identity(3,3) + (2*G0+G)*pow(dt_seg,2)/6.;
        MatrixXd PHI_rv = MatrixXd::Identity(3,3)*dt_seg +(G0+G)*pow(dt_seg,3)/12.;
        MatrixXd PHI_vr = (G0+G)*dt_seg/2.;
        MatrixXd PHI_vv = MatrixXd::Identity(3,3) + (2*G0+G)*pow(dt_seg,2)/6.;

        MatrixXd STM(6,6); STM << PHI_rr, PHI_rv, PHI_vr, PHI_vv;
        P1 = STM * _P0 * STM.transpose() + Qd_k;

        // update
        _r0 = r1_mean;
        _v0 = v1_mean;
        _P0 = P1;
    }

    //------------------------------------------//

    // Vector3d _r0 = r0_mean;
    // Vector3d _v0 = v0_mean;
    // MatrixXd _P0 =P0;
    
    // for (int k=0; k<nSub; k++)
    // {
    //     MatrixXd STM;
    //     propagateKEP_U(_r0, _v0, dt_sec/double(nSub), mu, r1_mean, v1_mean, STM);
    //     P1 = STM * _P0 * STM.transpose() + Qd_k;

    //     // cout << "[" << P1.block(0,0,2,2) << "]" << endl;

    //     // update
    //     _r0 = r1_mean;
    //     _v0 = v1_mean;
    //     _P0 = P1;

        
    // } 
    // // cout << "-----" << endl;

}                    