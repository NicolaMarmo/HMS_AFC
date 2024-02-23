#pragma once
#include <Eigen/Dense>

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;


 /* 6. UNIVERSAL VARIABLES FORMULATION */


//This function evaluates the Stumpff function S(z)
inline double stumpS(const double z)
{
    /*
     % This function evaluates the Stumpff function S(z) according
     % to Equation 3.49.
     %
     % z - input argument
     % s - value of S(z)
     %
     % User M - functions required : none
     % ------------------------------------------------------------*/
    double s;
    
    if (z > 0)
        s = (sqrt(z) - sin(sqrt(z))) / pow(sqrt(z), 3);
    else if (z < 0)
        s = (sinh(sqrt(-z)) - sqrt(-z)) / pow(sqrt(-z), 3);
    else
        s = 1. / 6.;
    
    return s;
}

//This function evaluates the Stumpff function C(z)
inline double stumpC(const double z)
{
    //    % This function evaluates the Stumpff function C(z) according
    //    % to Equation 3.50.
    //    %
    //    % z - input argument
    //    % c - value of C(z)
    //    %
    //    % User M - functions required : none
    //    % ------------------------------------------------------------
    double c;
    if (z > 0)
        c = (1 - cos(sqrt(z))) / z;
    else if (z < 0)
        c = (cosh(sqrt(-z)) - 1) / (-z);
    else
        c = 0.5;
    
    return c;
}

//Evaluates the Lagrange coefficients F,G,Fdot,Gdot
inline void LagrangeCoeff(const double r1, const double r2, const double p, const double mu, const double Danu, 
            double& F, double& G, double& Fdot, double& Gdot)
{
    //INPUT:
    //    r1: initial radius
    //    r2: final radius
    //    p: semi-latus rectum
    //    mu: gravitational constant
    //    Danu: true anomaly variation from r1 to r2
    
    F = 1. - r2/p*(1. - cos(Danu));
    G = r1*r2*sin(Danu)/(sqrt(p*mu));
    Fdot = sqrt(mu/p)*((1. - cos(Danu))/sin(Danu))*((1. - cos(Danu))/p - 1./r2 - 1./r1);
    Gdot = 1. - r1/p*(1 - cos(Danu));
}

//This function uses Newton's method to solve the universal Kepler equation for the universal anomaly.
inline double kepler_U(const double mu, const double dt, const double ro, const double vro, const double ua)
{
    // This function uses Newton's method to solve the universal
    // Kepler equation for the universal anomaly.
    //
    // mu - gravitational parameter(km^3 / s^2)
    // x - the universal anomaly(km^0.5)
    // dt - time since x = 0 (s)
    // ro - radial position(km) when x = 0
    // vro - radial velocity(km / s) when x = 0
    // ua - reciprocal of the semimajor axis(1 / km)
    // z - auxiliary variable(z = a*x^2)
    // C - value of Stumpff function C(z)
    // S - value of Stumpff function S(z)
    // n - number of iterations for convergence
    // nMax - maximum allowable number of iterations
    //
    // User M - functions required : stumpC, stumpS
    // ------------------------------------------------------------
    
    //%...Set an error tolerance and a limit on the number of
    //% iterations :
    const double error = 1.e-8;
    const int nMax = 1000;
    //%...Starting value for x:
    double x = sqrt(mu)*abs(ua)*dt;
    //%...Iterate on Equation 3.62 until convergence occurs within the error tolerance :
    int n = 0;
    double ratio = 1;
    
    while ((abs(ratio) > error) & (n <= nMax))
    {
        n = n + 1;
        double C = stumpC(ua*x*x);
        double S = stumpS(ua*x*x);
        double F = ro*vro / sqrt(mu)*x*x*C + (1. - ua*ro)*pow(x,3)*S + ro*x - sqrt(mu)*dt;
        double dFdx = ro*vro / sqrt(mu)*x*(1. - ua*x*x*S) + (1. - ua*ro)*x*x*C + ro;
        ratio = F / dFdx;
        x = x - ratio;
    }
    
    // if (n > nMax)
    // /*...Deliver a value for x, but report that nMax was reached */
    //     cout << "Number of iterations of Kepler''s equation = " << n << " > " << nMax << endl;
    
    
    return x;
}

//This function calculates the Lagrange f and g coefficients.
inline void lagrangefg_kepU(const double mu, const double x, const double t, const double ro, const double ua,
            double& f, double& g)
{
    //function[f, g] = f_and_g(x, t, r0, a)
    //%
    //% This function calculates the Lagrange f and g coefficients.
    //%
    //% mu - the gravitational parameter(km^3 / s^2)
    //% ua - reciprocal of the semimajor axis(1 / km)
    //% r0 - the radial position at time t(km)
    //% t - the time elapsed since t(s)
    //% x - the universal anomaly after time t(km^0.5)
    //% f - the Lagrange f coefficient(dimensionless)
    //% g - the Lagrange g coefficient(s)
    //%
    //% User M - functions required : stumpC, stumpS
    //% ------------------------------------------------------------
    
    double x2 = x*x;
    double z = ua*x2;
    //%...Equation 3.66a e  3.66b:
    f = 1 - x2 / ro*stumpC(z);
    g = t - 1. / sqrt(mu)*pow(x,3)*stumpS(z);
}

//This function calculates the Lagrange fdot and gdot coefficients.
inline void lagrangefDotgDot_kepU(const double mu, const double x, const double r, const double ro, const double ua,
            double& fdot, double& gdot)
{
    //function[fdot, gdot] = fDot_and_gDot(x, r, ro, a)
    //%
    //% This function calculates the time derivatives of the
    //% Lagrange f and g coefficients.
    //%
    //% mu - the gravitational parameter(km�3 / s�2)
    //% a - reciprocal of the semimajor axis(1 / km)
    //% ro - the radial position at time t(km)
    //% t - the time elapsed since initial state vector(s)
    //% r - the radial position after time t(km)
    //% x - the universal anomaly after time t(km�0.5)
    //% fDot - time derivative of the Lagrange f coefficient(1 / s)
    //% gDot - time derivative of the Lagrange g coefficient
    //% (dimensionless)
    //%
    //% User M - functions required : stumpC, stumpS
    //% ------------------------------------------------------------
    
    double z = ua*x*x;
    //%...Equation 3.66c: 3.66d :
    fdot = sqrt(mu) / r / ro*(z*stumpS(z) - 1)*x;
    gdot = 1. - x*x / r*stumpC(z);
}

//Propagation through universal elements. 
//This function computes the state vector (R, V) from the initial state vector (R0, V0) and the elapsed time dt_sec.
inline void propagateKEP_U(const Vector3d& R0, const Vector3d& V0, const double dt_sec, const double mu,
            Vector3d& R, Vector3d& V, MatrixXd& STM)
{
    //function[R, V] = rv_from_r0v0(R0, V0, t)
    //% This function computes the state vector(R, V) from the
    //% initial state vector(R0, V0) and the elapsed time.
    //%
    //% mu - gravitational parameter(km^3 / s^2)
    //% R0 - initial position vector(km)
    //% V0 - initial velocity vector(km / s)
    //% t - elapsed time(s)
    //% R - final position vector(km)
    //% V - final velocity vector(km / s)
    //%
    //% User M - functions required : kepler_U, f_and_g, fDot_and_gDot
    //% ------------------------------------------------------------
    
    //%...Magnitudes of R0 and V0 :
    double r0 = R0.norm();
    double v0 = V0.norm();
    //%...Initial radial velocity :
    double vr0 = R0.dot(V0) / r0;
    
    //%...Reciprocal of the semimajor axis(from the energy equation) :
    double alpha = 2 / r0 - v0*v0 / mu;
    
    //%...Compute the universal anomaly :
    double x = kepler_U(mu, dt_sec, r0, vr0, alpha);
    
    //%...Compute the f and g functions, derivatives of f and g
    double f, g;
    lagrangefg_kepU(mu, x, dt_sec, r0, alpha, f, g);
    
    //%...Compute the final position vector :
    R = f*R0 + g*V0;
    
    //%...Compute the magnitude of R :
    double r = R.norm();
    
    //%...Compute the derivatives of f and g
    double fdot, gdot;
    lagrangefDotgDot_kepU(mu, x, r, r0, alpha, fdot, gdot);
    
    //%...Compute the final velocity:
    V = fdot*R0 + gdot*V0;
    
    if (fabs(dt_sec) <= 1e-10)
    {
        R = R0;
        V = V0;
    }

    //% Compute the State Transition Matrix: [r;v] = STM * [r0;v0]
    STM = MatrixXd::Zero(6,6);
    STM << MatrixXd::Identity(3,3)*f, MatrixXd::Identity(3,3)*g,MatrixXd::Identity(3,3)*fdot, MatrixXd::Identity(3,3)*gdot;   
}

inline void propagateKEP_U(const Vector3d& R0, const Vector3d& V0, const double dt_sec, const double mu,
    Vector3d& R, Vector3d& V)
{
    MatrixXd STM;
    propagateKEP_U(R0, V0, dt_sec, mu, R, V, STM);
}