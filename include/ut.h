#pragma once
#include<eigen3/Eigen/Dense>

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

#include<iostream>
using namespace std;

#include "twoBP_lite.h"

void ut_propagation(int n);


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