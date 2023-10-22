#pragma once
#include<eigen3/Eigen/Dense>

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

#include<iostream>
using namespace std;

#include "twoBP_lite.h"
#include "C_sim_ode.h"
#include "ode.hpp"
#include "rk4_vector.h"

class C_UT
{
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





void propagate_UT(Vector3d r0_mean, Vector3d v0_mean, MatrixXd P0, double dt_sec, double mu, MatrixXd Qd_k, int Nsub,
                        Vector3d& r1_mean, Vector3d& v1_mean, MatrixXd& P1);

class EoM_Kepler{
    public:
    double mu;
    Vector3d T;
    double C;
    double tau;
    EoM_Kepler(double mu, Vector3d T, double C, double tau) : mu(mu), T(T), C(C), tau(tau) {};

    VectorXd operator()(double t, const VectorXd& x){
        VectorXd dxdt(7);
        Vector3d r = x.segment(0,3);
        double m = x(6);
        Vector3d a = -(mu*r)/(r.norm()*r.norm()*r.norm()) + T/m;
        double mp = -T.norm()/C;

        dxdt << tau*x.segment(3,3), tau*a, tau*mp;

        return dxdt;
    }
};

class EoM_CR3BP{
    public:
    double mu;
    double tau;
    EoM_CR3BP(double mu, double tau) : mu(mu), tau(tau) {};                   

    VectorXd operator()(double t, const VectorXd& X){
        VectorXd dxdt(6);
        double x = X(0); double y = X(1); double z = X(2); 
        double r1 = sqrt(pow((x + mu), 2) + pow(y, 2) + pow(z, 2)); double r1cube = pow(r1, 3);
        double r2 = sqrt(pow((x + mu - 1), 2) + pow(y, 2) + pow(z, 2)); double r2cube = pow(r2, 3);
        double vxp = 2*X(4) + x - (1 - mu)*(x + mu)/r1cube - mu*(x + mu - 1)/r2cube;
        double vyp = -2*X(3) + y - (1 - mu)*y/r1cube - mu*y/r2cube;
        double vzp = -(1 - mu)*z/r1cube - mu*z/r2cube;

        dxdt << tau*X.segment(3,3), tau*vxp, tau*vyp, tau*vzp;

        return dxdt;
    }
};

class EoM_CR3BP_ODE{
    public:
    double mu;
    double tau;
    EoM_CR3BP_ODE(double mu, double tau) : mu(mu), tau(tau) {};                   

    void operator()(double t, double X[], double dXdt[]){
        double x = X[0]; double y = X[1]; double z = X[2]; 
        double r1 = sqrt(pow((x + mu), 2) + pow(y, 2) + pow(z, 2)); double r1cube = pow(r1, 3);
        double r2 = sqrt(pow((x + mu - 1), 2) + pow(y, 2) + pow(z, 2)); double r2cube = pow(r2, 3);
        double vxp = 2*X[4] + x - (1 - mu)*(x + mu)/r1cube - mu*(x + mu - 1)/r2cube;
        double vyp = -2*X[3] + y - (1 - mu)*y/r1cube - mu*y/r2cube;
        double vzp = -(1 - mu)*z/r1cube - mu*z/r2cube;

        dXdt[0] = tau*X[3];
        dXdt[1] = tau*X[4];
        dXdt[2] = tau*X[5];
        dXdt[3] = tau*vxp;
        dXdt[4] = tau*vyp;
        dXdt[5] = tau*vzp;
    }
};