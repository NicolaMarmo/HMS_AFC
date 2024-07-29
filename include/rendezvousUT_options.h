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
    vector<Vector3d> r0_RV_vector, v0_RV_vector;
    vector<double> tfin_vector;

    vector<int> FB_Legs;
    
// Mission Switch
    string output_folder;
    string firstguess_folder;
    bool E_Pf_constraint;
    bool E_DV_cstr;

    bool v_RV_free;
    bool DV_RV_double;
    bool E_nav_std;
    bool Fixed_ToF_Leg;
    bool Limited_ToF;
    bool Equal_Segment_ToF;
    
    double mu_primary, mu_FB, SoI_R;
    double tfin1, tfin2, tfin3;
    double Max_ToF;

    double r_min;

    Vector3d r0, v0, rf, vf, rRV, vRV, r0_RV, v0_RV;

    double sigma2_r0, sigma2_v0, sigma2_rf, sigma2_vf, sigma2_rRV, sigma2_vRV;

    double DVtot_single_max;   

    int Qd_level;

    // Nondimensionalization
    double rconv;
    double vconv;
    double aconv;
    double tconv;

public:
    //costruttore
    Options_rendezvousUT_t(){};
    Options_rendezvousUT_t(const std::string &filename);
    
    //display
    void disp();
    // void print(const std::string &filename) const;
    void emit(const std::string &filename);

};


