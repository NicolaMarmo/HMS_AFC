#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <yaml-cpp/yaml.h>
#pragma warning(disable: 4996)
#pragma warning(disable: 4251)
#pragma warning(disable: 4275)
 
#include<eigen3/Eigen/Dense>
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
 
#include<vector>

using namespace std;

class Options_rendezvousUT_t
{
private:
    YAML::Node options_yaml;

public:
    
// Solver
    //int nSeg;
    int nSeg1, nSeg2, nSeg3;                    
    int nLeg;

    int sf, Max_iter;

    vector<int> nSeg_vector;
    vector<Vector3d> rRV_vector, vRV_vector, r0_RV_vector, v0_RV_vector;
    vector<double> tfin_vector;

    vector<int> FB_Legs;
    
// Mission Switch
    string output_folder;
    string firstguess_folder;
    bool E_Pf_constraint;
    bool E_DV_cstr;
    bool E_uniform_time;
    bool E_rendezvous;
    bool v_RV_free;
    bool DV_RV_double;
    // bool E_fakeUT;
    bool E_nav_std;
    bool objective_std;
    bool Fixed_ToF_Leg;
    bool Limited_ToF;
    int obj_func_switch;               // 0: energy, 1 = fuel;, 2 = min_trace_cov
    int DVstd_model;

    int cov_collocation_mode;     //0 = tutta, 1 = tiangolare (L), 2 = Simmetrica (S), 3 = Xcorr
    
    int cov_propagation_mode;
    int cov_propagation_nSub;
    
    double amu_dim, mu_FB;
    double tfin1, tfin2, tfin3;
    double Max_ToF;

    double r_min;

    Vector3d r0, v0, rf, vf, rRV, vRV, r0_RV, v0_RV;

    double sigma_r0;
    double sigma_v0;
    double sigma_rf_des;
    double sigma_vf_des;
    double sigma_rRV;
    double sigma_vRV;
    double DVtot_max;  
    double DVtot_single_max;   

    double Qd_r, Qd_v;
    int Qd_level;
    double nav_sigma_r;  //std^2
    double nav_sigma_v;  //std^2
 
    // bounds
    double px_lb, px_ub; 
    double py_lb, py_ub;
    double pz_lb, pz_ub;
    double vx_lb, vx_ub;
    double vy_lb, vy_ub;
    double vz_lb, vz_ub;

    // Nondimensionalization
    double rconv;
    double vconv;
    double aconv;
    double tconv;
    double amrif;

public:
    //costruttore
    Options_rendezvousUT_t(){};
    Options_rendezvousUT_t(const std::string &filename);
    
    //display
    void disp();
    // void print(const std::string &filename) const;
    void emit(const std::string &filename);

};


