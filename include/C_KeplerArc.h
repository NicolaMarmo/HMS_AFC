#pragma once
#include "ut.h"

#include<vector>
#include <numeric> //accumulate

#include<eigen3/Eigen/Dense>
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;



class C_KeplerArc
{
    public:
    Vector3d r0;
    Vector3d v0;
    double t0;
    double tf;

    bool final_state;
    Vector3d _rf, _vf;

    C_KeplerArc(const Vector3d& r0, const Vector3d& v0, const double& t0, const double&dt) : r0(r0), v0(v0), t0(t0)
    {
        tf = (t0 + dt);
        final_state = false;
    };

    void get_final_state(Vector3d& rf, Vector3d& vf);
    
    void print(ostream& os, int npt, bool header_on=false);


  
};


class C_KeplerMultiArc
{
    public:

    vector<C_KeplerArc> v_keplerArcs;


    C_KeplerMultiArc(){};
    /** 
     * 
     * Argument size:
     *      v_t[N+1]
     *      MatrixXd DVs[3, N+1]
     */
    C_KeplerMultiArc(Vector3d _r0, Vector3d _v0, VectorXd v_t, MatrixXd DVs);

    void add_keplerArc(C_KeplerArc arc);
    void print(ostream& os, int npt);

};



