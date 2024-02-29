#pragma once
#include <Eigen/Dense>

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

#include <iostream>
using namespace std;

#include "twoBP_lite.h"
#include "rk4_vector.h"
#include "ode.hpp"
#include "C_sim_ode.h"

void ut_propagation(int n);


class C_UT{
    public:

    int n;
    VectorXd W;

    C_UT(int n, double w0=0.) : n(n)
    {

        W = VectorXd(1+2*n);
        
        W(0) = w0;
        for (int j=1; j<1+2*n; j++)
            W(j) = (1.-w0)/(2*n);

        //DEBUG
        // cout << "W = " << W.transpose() << endl;

    };

    void eval_sigmaPoints(VectorXd xMean, MatrixXd Pxx, MatrixXd& XX)
    {
        int n = xMean.size();

        Eigen::MatrixXd L( (Pxx*n/(1.-W(0))).llt().matrixL() );

        XX = MatrixXd(n,1+2*n);
        XX.col(0)=xMean;
        for (int j=0; j<n; j++)
            XX.col(1+j) = xMean + L.col(j);

        for (int j=0; j<n; j++)
            XX.col(1+n+j) = xMean - L.col(j);
            
    }


    void inverse_sigmaPoints(MatrixXd YY, VectorXd& y_mean, MatrixXd& Pyy)
    {
           
        // mean
        // cout << "n cols of YY: " << YY.cols() << endl;
        // cout << YY.col(0) << endl;
        //   cout << W(0) << endl;

        y_mean = VectorXd::Zero(n);
        for (int j=0; j<2*n+1; j++)
            y_mean = y_mean + W[j]*YY.col(j);

        // covariance
        Pyy = MatrixXd::Zero(n,n);
        for (int j=0; j<2*n+1; j++)
            Pyy += W(j)*(YY.col(j)-y_mean)*(YY.col(j)-y_mean).transpose();

    }

};

void propagate_kepler_UT(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k,
                        Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1);

void propagate_kepler_UT_MS(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k,
                        Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1);
/**
 * @brief Propaga usando la state transition matrix
 * 
 */
void propagate_kepler_STM(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k,
                         Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1, int nSub=1);

void Pf_STM(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k, 
            MatrixXd& P1);                         

// class EoM_Kepler{
//     public:
//     double mu;
//     double tau;
//     EoM_Kepler_STM(double mu, double tau): mu(mu), tau(tau) {};

//     VectorXd operator()(double t, const VectorXd& x){
//         VectorXd dxdt(6);
//         Vector3d r = x.segment(0,3);
//         Vector3d a = -(mu*r)/(r.norm()*r.norm()*r.norm());

//         dxdt << tau*x.segment(3,3), tau*a;

//         return dxdt;
//     }
// };

class EoM_Kepler_STM{
    public:
    double mu;
    double tau;

    EoM_Kepler_STM(double mu, double tau) : mu(mu), tau(tau) {};

    VectorXd operator()(double t, const VectorXd& x){
        VectorXd dxdt(42);
        Vector3d r = x.segment(0, 3);
        Vector3d v = x.segment(3, 3);

        Vector3d a = -(mu*r)/(r.norm()*r.norm()*r.norm());

        // Compute the STM components
        MatrixXd A(6, 6);
        A.setZero();
        A.block<3, 3>(0, 3) = Matrix3d::Identity();
        A.block<3, 3>(3, 0) = (3*mu*r*r.transpose() - mu*pow(r.norm(), 2)*Matrix3d::Identity())/pow(r.norm(), 5);

        // Create the state transition matrix Phi
        MatrixXd Phi(6, 6), Phip(6,6);
        Phi = x.tail(36).reshaped(6, 6);
        Phip = A*Phi*tau;

        // Combine state derivatives and STM
        dxdt.segment(0, 6) << tau*v, tau*a;
        dxdt.segment(6, 36) = Phip.reshaped();

        return dxdt;
    }
};

class EoM_Kepler_ODE{
    public:
    double mu;
    double tau;
    EoM_Kepler_ODE(double mu, double tau) : mu(mu), tau(tau) {};                   

    void operator()(double t, double X[], double dXdt[]){
        double rcube = pow((pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2)), 1.5); 

        dXdt[0] = tau*X[3];
        dXdt[1] = tau*X[4];
        dXdt[2] = tau*X[5];
        dXdt[3] = -tau*(mu*X[0])/rcube;
        dXdt[4] = -tau*(mu*X[1])/rcube;
        dXdt[5] = -tau*(mu*X[2])/rcube;
    }
};

class EoM_Kepler_ODE_STM{
    public:
    double mu;
    double tau;
    EoM_Kepler_ODE_STM(double mu, double tau) : mu(mu), tau(tau) {};                   

    void operator()(double t, double X[], double dXdt[]){
        Vector3d r; r << X[0], X[1], X[2];
        double r3 = pow(r.norm(), 3); 
        double r2 = pow(r(0), 2) + pow(r(1), 2) + pow(r(2), 2); 
        MatrixXd Phi(6, 6), Phip(6, 6);
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                Phi(i, j) = X[6*i + j + 6];
            }
        }

        MatrixXd A(6, 6);
        A.setZero();
        A.block<3, 3>(0, 3) = Matrix3d::Identity();
        A.block<3, 3>(3, 0) = mu*(3*r*r.transpose() - r2*MatrixXd::Identity(3, 3))/pow(r.norm(), 5);

        Phip = A*Phi*tau;

        dXdt[0] =  tau*X[3];
        dXdt[1] =  tau*X[4];
        dXdt[2] =  tau*X[5];
        dXdt[3] = -tau*(mu*X[0])/r3;
        dXdt[4] = -tau*(mu*X[1])/r3;
        dXdt[5] = -tau*(mu*X[2])/r3;
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                dXdt[6*i + j + 6] = Phip(i, j);
            }
        }
    }
};

void Propagate_P_FB(const MatrixXd& P0, const Vector3d& v_infm, const Vector3d& v_infp, const double SoI_R, const double mu,
    MatrixXd Pf);

void ToF_Hyperbola(const Vector3d r_infm, const Vector3d v_infm, const double mu, double& DT);