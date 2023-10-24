#pragma once

#include "ut.h"
#include "C_KeplerArc.h"
#include "rendezvousUT_options.h"
#include "C_sim_ode.h"
#include "ode.hpp"

#include "aux_covariance.h"

#include<worhp/worhp.h>
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

void test_RR_HMS();

class C_prb_RR_HMS
{
    public:
    double rconv, vconv, tconv, aconv, amrif;
    double prconv; 

    const double muSun = 132712440018;      //km^3/s^2
    const double muEarth = 398600.4;        //km^3/s^2
    const double muMars = 42828;            //km^3/s^2
    const double REarth = 6378.136;         //km
    const double RMars = 3389.5;            //km

    double t0;
    double tf1;
    double tf2;

    MatrixXd P0;
    Vector3d r0, v0; 
    Vector3d rf_des, vf_des;
    Vector3d rRV, vRV;
    Vector3d r_G, v_G;
    Vector3d r0_X, v0_X; 

    double sigma_r_G, sigma_v_G;

    double mass0;
    double THR, C;
    MatrixXd Qd_k;
    MatrixXd RR_k;  // covariance of navigation error

    ifstream opt_Qd;
    string line;
    string el;
    VectorXd X;

    //
    bool E_covariance_propagation;   //not used
    double kstd; // coefficiente che tiene conto del "DV stocastico"        DV = DV_det + kstd * DV_std

    // opzioni
    int nSeg;
    int nSeg1;       // n. di segmenti in cui è divisa la traiettoria
    int nSeg2;
    int nLeg;
    int sf;
    const int Pdim = 6;

    //double *tf_vector_p;
    vector<double> tf_vector;
    vector<int> nSeg_vector;
    vector<Vector3d> rRV_vector;
    vector<Vector3d> vRV_vector;
    
    const int nRes = 6;     //n. di equazioni di continuità
    const int nVars = 1 + 6 + 3 + 18; // YP [nVars x nSeg]    
    Options_rendezvousUT_t opz;


    /*
	 * WORHP data structures
     *
     * OptVar contains the variables, constraint values and multipliers.
     * Workspace encapsulates all workspace and working variables for WORHP.
     * Params contains all WORHP parameters.
     * Control contains things for reverse communication flow control.
     */
    OptVar    opt;
    Workspace wsp;
    Params    par;
    Control   cnt;


    C_prb_RR_HMS(string fname_opz);
    int setup_worph(double *Xguess);
    void evalFG(double *X, double& F, double *G, double ScaleObj);


    void custom_iteration_process(double *X_current);


    /*
     * WORHP Reverse Communication loop.
     * In every iteration poll GetUserAction for the requested action, i.e. one
     * of {callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif}.
     *
     * Make sure to reset the requested user action afterwards by calling
     * DoneUserAction, except for 'callWorhp' and 'fidif'.
     */
    int solve()
    {
        double F_last;
        double *G_last;
        G_last = new double[opt.m];

        // opt.newX=true;
        while (cnt.status < TerminateSuccess && cnt.status > TerminateError)
        {
            /*
         * WORHP's main routine.
         * Do not manually reset callWorhp, this is only done by the FD routines.
         */
            if (GetUserAction(&cnt, callWorhp))
            {
                Worhp(&opt, &wsp, &par, &cnt);
                // No DoneUserAction!
            }

            /*
         * Show iteration output.
         * The call to IterationOutput() may be replaced by user-defined code.
         */
            if (GetUserAction(&cnt, iterOutput))
            {
                IterationOutput(&opt, &wsp, &par, &cnt);

                custom_iteration_process(opt.X);


                DoneUserAction(&cnt, iterOutput);
            }

            /*
         * Evaluate the objective function.
         * The call to UserF may be replaced by user-defined code.
         */
            if (GetUserAction(&cnt, evalF))
            {
                evalFG(opt.X, opt.F, opt.G, wsp.ScaleObj);
                DoneUserAction(&cnt, evalF);
            }

            /*
         * Evaluate the constraints.
         * The call to UserG may be replaced by user-defined code.
         */
            if (GetUserAction(&cnt, evalG))
            {
                // // UserG(&opt, &wsp, &par, &cnt);
                // if (opt.newX)
                // {
                //     evalFG(opt.X, F_last, G_last);
                // }
                // // opt.F = F_last;
                // opt.G = G_last;
                evalFG(opt.X, opt.F, opt.G, wsp.ScaleObj);
                DoneUserAction(&cnt, evalG);
            }

            /*
         * Use finite differences with RC to determine derivatives
         * Do not reset fidif, this is done by the FD routine.
         */
            if (GetUserAction(&cnt, fidif))
            {
                WorhpFidif(&opt, &wsp, &par, &cnt);
                // No DoneUserAction!
            }


        }

        // Translate the WORHP status flag into a meaningful message.
        StatusMsg(&opt, &wsp, &par, &cnt); 
        return EXIT_SUCCESS;
    }

    void clean()
    {
        // Deallocate all data structures.
        // Data structures must not be accessed after this call.
        WorhpFree(&opt, &wsp, &par, &cnt);
    }

    
    /**
     * @brief inizializzazione piana,
     * 
     * @param Xguess [nVars*(nSeg+1)]
     */
    void init_lin(double *Xguess);    


    void save_sol(string fname="../results/fullsol.dat", bool skipAsking=false)
    {       
        string ans = "S";
        if (not(skipAsking)){
            //cout << "Salvare la soluzione?" << endl;
            //cin >> ans;
        }
        if (ans.compare("S") == 0 || ans.compare("s") == 0)
        {
            ofstream out_fullsol(fname);
            for (int iLeg=0; iLeg < nLeg; iLeg++){
                nSeg = nSeg_vector[iLeg];
                for (int k=0; k<nSeg+1; k++) // Warning!
                {
                    int i0 = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + k*nVars;
                    //int i0 = iLeg*(nSeg+1)*nVars + k*nVars;
                    for (int j = 0; j < nVars; j++)
                        out_fullsol << fixed << setw(14) << setprecision(8) << opt.X[i0 + j] << " ";
                    out_fullsol << endl;
                }
                //out_fullsol << endl << endl;
            }
            out_fullsol.close();
        }
    };

    void load_sol(string fname, double* Xguess)
    {
        std::string ans;
        cout << "caricare il file soluzione? (S/n)" << endl;
        cin >> ans;
        if (ans.compare("S") == 0 || ans.compare("s") == 0)
        {
            vector<double> v_xGuess;
            std::ifstream infile(fname);
            std::string line;
            while (std::getline(infile, line))
            {
                std::istringstream iss(line);
                double num;
                while (iss >> num) 
                    v_xGuess.push_back(num);

                // process pair (a,b)
            }
            
            for (int i=0; i<v_xGuess.size(); i++){
                Xguess[i] = v_xGuess[i];
            }
        }
    }

};





