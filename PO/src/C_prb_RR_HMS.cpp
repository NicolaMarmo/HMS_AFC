#include "../include/C_prb_RR_HMS.h"
#include<string>
#include<sstream>
#include <ctime>

 void C_prb_RR_HMS::init_lin(double *Xguess){
    int krev = 0;

    // Default init: X=0
    for (int ii = 0; ii < nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); ii++) 
        Xguess[ii] = 0.;

    for (int iLeg = 0; iLeg<nLeg; iLeg++){
        int nSeg = nSeg_vector[iLeg];
        double tf = tfin_vector[iLeg];
        double dt = tf/nSeg;
        VectorXd tspan = VectorXd::LinSpaced(nSeg + 1, t0, tfin_vector[iLeg]);
        
        MatrixXd stateGuess(6, nSeg + 1);

        VectorXd s0(6), sf(6);
        if (iLeg == 0){
            s0 << r0, v0;
        }
        else{
            s0 << rRV_vector[iLeg - 1], vRV_vector[iLeg - 1];
        }

        for (int k = 0; k < nSeg + 1; k++){
            C_simulazione_ODE sim(6, 50, 1, 1e-8, 1e-10);
            sim.SetY0_Sim(s0);
            EoM_CR3BP_ODE EoM(mu, k*dt);
            sim.Start_Sim(EoM);
            sf = sim.Y.rightCols(1);

            stateGuess(0, k) = sf(0);
            stateGuess(1, k) = sf(1);
            stateGuess(2, k) = sf(2);
            stateGuess(3, k) = sf(3);
            stateGuess(4, k) = sf(4);
            stateGuess(5, k) = sf(5);
        }

        int iX = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg);
        for (int k=0; k<nSeg+1; k++){
            int i0 = nVars*k;

            // state
            for (int j=0; j<6; j++)
                Xguess[iX + i0 + j] = stateGuess(j, k);
                        
            // DV [null]            

            // KK [null]
        }
    }
};

/**
 * @brief Construct a new problem:
 *          collocazione dello stato, propagazione della covarianza
 *      
 * 
 * @param fname_opz 
 */

C_prb_RR_HMS::C_prb_RR_HMS(string fname_opz){
    //--------------------------------------------
    // Carica le opzioni
    //
    opz = Options_rendezvousUT_t(fname_opz);
    
    rconv = opz.rconv;
    vconv = opz.vconv;
    tconv = opz.tconv;
    aconv = opz.aconv;
    
    nLeg = opz.nLeg;
    sf = opz.sf;

    Nsub = opz.Nsub;

    mu = opz.mu;

    cout << "rconv = " << rconv << endl;
    cout << "vconv = " << vconv << endl;
    cout << "tconv = " << tconv << endl;
    cout << "aconv = " << aconv << endl;

    t0 = 0;

    nSeg_vector = opz.nSeg_vector;
    rRV_vector = opz.rRV_vector;
    vRV_vector = opz.vRV_vector;
    tfin_vector = opz.tfin_vector;

    r0 = opz.r0;
    v0 = opz.v0;
    rf = opz.rf;
    vf = opz.vf;

    P0 = MatrixXd::Zero(6, 6);
    P0.diagonal() << opz.sigma_r0, opz.sigma_r0, opz.sigma_r0, opz.sigma_v0, opz.sigma_v0, opz.sigma_v0;
    
    Qd_k = MatrixXd::Zero(Pdim, Pdim); 
    if (opz.planar)
        Qd_k.diagonal() << opz.Qd_r, opz.Qd_r, 0., opz.Qd_v, opz.Qd_v, 0;
    else {
        //Qd_k.diagonal() << opz.Qd_r, opz.Qd_r, opz.Qd_r, opz.Qd_v, opz.Qd_v, opz.Qd_v;
        //Qd_k.diagonal() << opz.Qd_r, opz.Qd_r, 0., opz.Qd_v, opz.Qd_v, 0;
        if (opz.Qd_level == 0)
            Qd_k = MatrixXd::Zero(Pdim, Pdim); 
        else{
            X = VectorXd::Zero(Pdim);
            string Qd_level_array[5] = {"Qds", "Qdm", "Qdl", "Qdxl", "Qdxxl"};
            opt_Qd.open ("../src/" + Qd_level_array[opz.Qd_level - 1] + ".dat", ios::in); 
            for (int i = 0; i < 6; i++){
                getline(opt_Qd, line);
                stringstream ss(line); 
                for (int j = 0; j < 6; j++){
                    ss >> el;
                    X(j) = stod(el);
                }
                Qd_k.row(i) = X;
            }
            opt_Qd.close();
            }
    }
    
    RR_k = MatrixXd::Zero(Pdim, Pdim); 
    if (opz.planar)
        RR_k.diagonal() << opz.nav_sigma_r, opz.nav_sigma_r, 0., opz.nav_sigma_v, opz.nav_sigma_v, 0;
    else
        RR_k.diagonal() << opz.nav_sigma_r, opz.nav_sigma_r, opz.nav_sigma_r, opz.nav_sigma_v, opz.nav_sigma_v, opz.nav_sigma_v;

    // const double beta = 0.05; // Pr{DV > DVmax} < 1 - beta
    kstd = sqrt(2*log(1./opz.beta) + sqrt(3.)); // circa 4.18 @ beta=0.05

    cout << "kstd:" << kstd << endl;

};

int C_prb_RR_HMS::setup_worph(double *Xguess){
    // Setup worph

    // Check Version of library and header files
    CHECK_WORHP_VERSION

    // Properly zeros everything, or else the following routines could get confused
    WorhpPreInit(&opt, &wsp, &par, &cnt);

    // Uncomment this to get more info on data structures
    //WorhpDiag(&opt, &wsp, &par, &cnt);

    /*
    * Parameter initialisation routine that must be called
    * when using ReadParamsNoInit instead of ReadParams.
    */
    int status;
    InitParams(&status, &par);

    /*
    * We can now set parameters that may be overruled by those in the
    * parameter file. This is useful for setting a non-default standard
    * parameter value that may still be overwritten.
    */
    par.NLPprint = 1; // Let's prefer the slim output format
                      // unless the parameter file says differently

    /*
    * Parameter XML import routine that does not reset
    * all parameters to default values (InitParams does this)
    */
    ReadParamsNoInit(&status, "../worhp.xml", &par);
    if (status == DataError || status == InitError)
    {
        return EXIT_FAILURE;
    }

    /*
    * WORHP data structure initialisation routine.
    * Calling this routine prior to WORHP is mandatory.
    * Before calling WorhpInit, set the problem and matrix dimensions as
    *
    * opt.n      = number of variables,
    * opt.m      = number of constraints (lin + nonlin, excluding box con's),
    * wsp.DF.nnz = nonzero entries of the objective function gradient,
    * wsp.DG.nnz = nonzero entries of the constraint Jacobian,
    * wsp.HM.nnz = nonzero entries of the Lagrange Hessian.
    *
    * Set nnz to 'WorhpMatrix_Init_Dense' to have WorhpInit allocate and
    * create a dense matrix structure appropriate for the matrix kind and
    * its dimensions. Setting it to its dense dimension achieves the same.
    */
    opt.n = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); // {r,vm,Pm,DVx, DVy, DVz, K} x (N+1)
    
    opt.m = 6*nLeg;              // vincoli all'arrivo
    opt.m += nRes*accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.);     //continuità r,v punti interni + covarianza    
    if (opz.E_Pf_constraint) opt.m += 6*nLeg;

    if (opz.E_DV_cstr) opt.m += accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg; //ciascun DV è vincolato

    // All derivatives for this problem have a sparse structure, so  set the amount of nonzeros here
    wsp.DF.nnz = WorhpMatrix_Init_Dense;  
    wsp.DG.nnz = WorhpMatrix_Init_Dense;  
    wsp.HM.nnz = WorhpMatrix_Init_Dense;  

    WorhpInit(&opt, &wsp, &par, &cnt);
    if (cnt.status != FirstCall)
    {
        std::cout << "Main: Initialisation failed." << std::endl;
        return EXIT_FAILURE;
    }

    /*
    * These pointers give access to the essential user data:
    *
    * opt.X[0] to opt.X[opt.n - 1]           : Optimisation variables
    * opt.Lambda[0] to opt.Lambda[opt.n - 1] : Multipliers for the constraints
    *                                          on X ("box constraints")
    * opt.G[0] to opt.G[opt.m - 1]           : Linear and nonlinear constraints
    * opt.Mu[0] to opt.Mu[opt.m - 1]         : Multipliers for the constraints on G
    *
    * Set initial values of X, Lambda and Mu here.
    * G need not be initialised.
    */
    for (int i = 0; i < opt.n; i++) opt.X[i] = Xguess[i];

    /*
    * Set lower and upper bounds on the variables and constraints.
    * Use +/-par.Infty to signal "unbounded".
    *
    * XL and XU are lower and upper bounds ("box constraints") on X.
    * GL and GU are lower and upper bounds on G.
    */

    // bound generici 
    for (int iLeg = 0; iLeg < nLeg; iLeg++){
        int nSeg = nSeg_vector[iLeg];
        for (int iSeg = 0; iSeg < nSeg + 1; iSeg++){
            int i0 = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + iSeg*nVars;
            //int i0 = iLeg*(nSeg+1)*nVars + iSeg*nVars;
            for (int j=0; j<3; j++){
                opt.XL[i0+j] = -2;
                opt.XU[i0+j] = +2;
            }
            for (int j=0; j<3; j++){
                opt.XL[i0+3+j] = -2;
                opt.XU[i0+3+j] = +2;
            }
    
            for (int j=0; j<3; j++){
                if (opz.E_DV_cstr){   
                    double LDV;
                    double UDV;
                    if ((opz.DV_RV_double == false) && iSeg == nSeg){
                        LDV = UDV = 0;
                    }
                    else{
                        LDV = -opz.DVtot_single_max;
                        UDV = +opz.DVtot_single_max;
                    }
                    // opt.XL[i0+6+j] = -opz.DVtot_max/double(nSeg+1);
                    // opt.XU[i0+6+j] = +opz.DVtot_max/double(nSeg+1);
                    opt.XL[i0+6+j] = LDV;
                    opt.XU[i0+6+j] = UDV;
                }
                else{
                    opt.XL[i0+6+j] = -1;
                    opt.XU[i0+6+j] = +1;
                }
            }

            // KK
            if (opz.E_Pf_constraint){
                double LKK;
                double UKK;
                if ((opz.DV_RV_double == false) && iSeg == nSeg){
                    LKK = UKK = 0;
                }
                else{
                    LKK = -200;
                    UKK = +200;
                }
                //if(iLeg > 0) LKK = UKK = 0;
                for (int iRow=0; iRow<3; iRow++){       
                    for (int iCol=0; iCol<6; iCol++){
                        opt.XL[i0+6+3+iRow*6+iCol] = LKK;
                        opt.XU[i0+6+3+iRow*6+iCol] = UKK;
                    }
                }
            }
            else{
                for (int iCol=0; iCol<3; iCol++){       
                    for (int iRow=0; iRow<6; iRow++){
                        opt.XL[i0+6+3+iCol*3+iRow] = opt.XU[i0+6+3+iCol*3+iRow] = 0;
                    }
                }
            }
        }

        // Condizioni iniziali assegnate
        if (iLeg == 0){
            r0_X = r0;
            v0_X = v0;
        }
        else{
            r0_X = rRV_vector[iLeg - 1];
            v0_X = vRV_vector[iLeg - 1];
        }

        int i0 = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg);
        for (int j=0; j<3; j++){
            opt.XL[i0 + j] = opt.XU[i0 + j] = r0_X(j);
        }
        if((opz.v_RV_free == false)||(iLeg == 0)){                                 
            for (int j=0; j<3; j++){
                opt.XL[i0 + 3+j] = opt.XU[i0 + 3+j] = v0_X(j);
            }  
        }

        // Vincoli: Residui + Condizioni al contorno
        // residuo r & v
        for (int iSeg = 0; iSeg < nSeg; iSeg++){
            int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + (accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + iSeg*(nRes);
            //int iG0 = iLeg*(nSeg*nRes + 6 + nSeg + 1) + iSeg*(nRes);
            if (opz.E_Pf_constraint) iG0 += 6*iLeg;

            //residuo_rv
            for (int j = 0; j < 6; j++){  
                opt.GL[iG0 + j] = opt.GU[iG0 + j] = 0; 
            }
        }

        // condizioni contorno finale
        int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + (accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + nSeg*(nRes);
        //int iG0 = iLeg*(nSeg*nRes + 6 + nSeg + 1) + nSeg*(nRes);
        if (opz.E_Pf_constraint) iG0 += 6*iLeg;

        for (int i = 0; i < 6; i++){
            opt.GL[iG0 + i] = opt.GU[iG0 + i] = 0; 
        }

        if (opz.obj_func_switch == 1){ 
            // covarianza di posizione finale: Pr - Prdes < 0
            if (opz.E_Pf_constraint){
                for (int j=0; j<6; j++)
                {
                    opt.GL[iG0 + 6 + j] = -par.Infty;   
                    opt.GU[iG0 + 6 + j] = 0;
                }
            } 
            // DVcstr: DV[k] - DV_max < 0
            if (opz.E_DV_cstr){
                int nPf_cstr = 0;
                if (opz.E_Pf_constraint) nPf_cstr = 6;
                for (int j = 0; j < nSeg + 1; j++){
                    opt.GL[iG0 + 6 + nPf_cstr + j] = -par.Infty; 
                    opt.GU[iG0 + 6 + nPf_cstr + j] = 0;
                }
            } 
        }
        else{
            cout << " opz.obj_func_switch = " << opz.obj_func_switch << " non definito" << endl;
        } 
    }
    par.FGtogether = true;
    par.UserDF = false;
    par.UserDG = false;
    par.UserHM = false;
    return 0;
}

void C_prb_RR_HMS::custom_iteration_process(double *X_current){
    return;
}

void C_prb_RR_HMS::evalFG(double *X, double &F, double *G, double ScaleObj){
    MatrixXd eye6x6 = MatrixXd::Identity(6, 6);
    C_UT ut(6);
    mu = opz.mu;
    MatrixXd P_f;
    MatrixXd P_k1m;
    Vector3d r_km;
    Vector3d v_km;                                     
    MatrixXd P_km;
    Vector3d v_kp;
    MatrixXd P_kp;
    Vector3d v_k1p;
    MatrixXd covDV;
    Vector3d DV_k;
    MatrixXd KK_k;
    MatrixXd KKaug_k(6, 6);  

    Vector3d r_k1_hat, v_k1m_hat;  MatrixXd P_k1m_hat;

    VectorXd Pf_ev_r(3);
    VectorXd Pf_ev_v(3);
    F = 0;
    
    for (int iLeg = 0; iLeg < nLeg; iLeg++){
        int nSeg = nSeg_vector[iLeg];
        vector<double> v_DVnorm(nSeg + 1); 
        vector<double> v_DVstd(nSeg + 1);

        double dt = (tfin_vector[iLeg] - t0)/nSeg;
    
        int iX = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg);
        //int iX = iLeg*(nSeg+1)*nVars;

        r_km = Eigen::Map<Vector3d>(X + iX, 3);  // Initial state is assigned by setting upper and lower boundaires in X equal to \tilde{x0} 
        v_km = Eigen::Map<Vector3d>(X + iX + 3, 3);    
        if (iLeg == 0){
            P_km = P0;  
        }        
        else{
            P_km = P_kp;
            if (opz.v_RV_free) v_km = v_kp;
        }                          
        DV_k = Eigen::Map<Vector3d>(X + iX + 6, 3);      
        KK_k = Eigen::Map<MatrixXd>(X + iX + 6 + 3, 3, 6); 
        KKaug_k << MatrixXd::Zero(3, 6), KK_k; 

        // apply DV
        v_kp = v_km + DV_k;
        // cout << KK_k << endl << endl;
        P_kp = (eye6x6 + KKaug_k) * P_km * (eye6x6 + KKaug_k).transpose();
        covDV = KK_k*P_km*KK_k.transpose();

        v_DVnorm[0] = DV_k.norm();  
        if(opz.DVstd_model == 1){
            v_DVstd[0] = sqrt(covDV.eigenvalues().real().maxCoeff());      // max autovalore
        } 
        else if(opz.DVstd_model == 0){
            v_DVstd[0] = sqrt(covDV.trace());
        }

        for (int iSeg = 0; iSeg < nSeg; iSeg++){   
            propagate_UT(r_km, v_kp, P_kp, dt, mu, Qd_k, Nsub, r_k1_hat, v_k1m_hat, P_k1m_hat);

            int i0 = iX + (iSeg+1)*nVars;
            Eigen::Map<Vector3d> r_k1m(X+i0, 3);  
            Eigen::Map<Vector3d> v_k1m(X+i0+3, 3);                                       
            P_k1m = P_k1m_hat;                                  // Hybrid: Covarianza propagata
            Eigen::Map<Vector3d> DV_k1(X+i0+6, 3);        
            Eigen::Map<MatrixXd> KK_k1(X+i0+6+3, 3, 6);        
            
            int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + (accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + iSeg*(nRes);
            //int iG0 = iLeg*(nRes*nSeg + 6 + nSeg + 1) + nRes*(iSeg); 
            if (opz.E_Pf_constraint) iG0 += 6*iLeg;

            for (int j=0; j<3; j++)
                G[iG0 + j] = r_k1m[j] - r_k1_hat[j];
            for (int j=0; j<3; j++)
                G[iG0 + 3 + j] = v_k1m[j] - v_k1m_hat[j];

            // apply DV
            v_k1p = v_k1m + DV_k1;
            MatrixXd KKaug_k1(6, 6); KKaug_k1 << MatrixXd::Zero(3, 6), KK_k1;  
            MatrixXd P_k1p = (eye6x6 + KKaug_k1) * P_k1m_hat * (eye6x6 + KKaug_k1).transpose();  
            MatrixXd cov_DV_k1 = KK_k1*P_k1m*KK_k1.transpose();    
            
            v_DVnorm[iSeg + 1] = DV_k1.norm();  
            if(opz.DVstd_model == 1){
                v_DVstd[iSeg + 1] = sqrt(cov_DV_k1.eigenvalues().real().maxCoeff());      // max autovalore
            } 
            else if(opz.DVstd_model == 0){
                v_DVstd[iSeg + 1] = sqrt(cov_DV_k1.trace());
            }
            // update
            r_km = r_k1m; // note that: r_km = r_kp
            v_kp = v_k1p;
            P_kp = P_k1p;
        }

        // obj function
        double DVtot = accumulate(v_DVnorm.begin(), v_DVnorm.end(), 0.);
        double DVstd_tot = accumulate(v_DVstd.begin(), v_DVstd.end(), 0.);
        //     F = DVtot;
        F += DVtot + DVstd_tot;   //* kstd

        // constraints
        //
        
        if (iLeg == nLeg - 1){
            r_G = opz.rf;
            v_G = opz.vf;
        }
        else{
            r_G = rRV_vector[iLeg];
            v_G = vRV_vector[iLeg];
        }
        
        // 1.1) Stato medio al tempo finale
        int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + (accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + nSeg*(nRes);
        //int iG0 = iLeg*(nRes*nSeg + 6 + nSeg + 1) + nRes*nSeg;
        if (opz.E_Pf_constraint) iG0 += 6*iLeg;

        if ((opz.v_RV_free)&&(iLeg < nLeg - 1)){
            for (int j=0; j<3; j++)
                G[iG0 + j] = r_km[j] - r_G[j];

            for (int j=0; j<3; j++)
                G[iG0 + 3 +j] = 0;
        }
        else{
            for (int j=0; j<3; j++)
                G[iG0 + j] = r_km[j] - r_G[j];

            for (int j=0; j<3; j++)
                G[iG0 + 3 +j] = v_kp[j] - v_G[j];
        }

        // 1.2) Covarianza al tempo finale
        if (opz.E_Pf_constraint){
            // Pf_ev_r = P_kp.block(0,0,3,3).eigenvalues().real();
            // Pf_ev_v = P_kp.block(3,3,3,3).eigenvalues().real();

            Pf_ev_r = P_k1m_hat.block(0,0,3,3).eigenvalues().real();
            Pf_ev_v = P_k1m_hat.block(3,3,3,3).eigenvalues().real();
            
            if (iLeg == nLeg - 1){
                sigma_r_G = opz.sigma_rf_des;
                sigma_v_G = opz.sigma_vf_des;
                }
            else{
                sigma_r_G = opz.sigma_rRV;
                sigma_v_G = opz.sigma_vRV;
            }
             
            for (int j=0; j<3; j++){
                //G[iG0 + 6 + j] = P_kp(j,j) - sigma_r_G; 
                if(opz.sigma_constraint == 1){
                    G[iG0 + 6 + j] = opz.sf*(P_k1m_hat(j,j) - sigma_r_G);
                }
                else{
                    G[iG0 + 6 + j] = opz.sf*(Pf_ev_r(j) - sigma_r_G);
                }
            }
            for (int j=3; j<6; j++){
                //G[iG0 + 6 + j] = P_kp(j,j) - sigma_v_G; 
                if(opz.sigma_constraint == 1){
                    G[iG0 + 6 + j] = opz.sf*(P_k1m_hat(j,j) - sigma_v_G);  
                }
                else{
                    G[iG0 + 6 + j] = opz.sf*(Pf_ev_v(j-3) - sigma_v_G); 
                }
            }
        }

        // 1.3) DV_cstr
        if (opz.E_DV_cstr){
            int nPf_cstr = 0;
            if (opz.E_Pf_constraint) nPf_cstr = 6;

            for (int k=0; k<nSeg+1; k++){
                // G[iG0+6+nPf_cstr + k] = v_DVnorm[k] + kstd*v_DVstd[k] - 1*opz.DVtot_max/(nSeg+1); 
                G[iG0+6+nPf_cstr + k] = v_DVnorm[k] + kstd*v_DVstd[k] - 1*opz.DVtot_single_max; 
            }
        }
    }
    F *= ScaleObj;
}
 
void test_RR_HMS(){
    MatrixXd eye6 = MatrixXd::Identity(6, 6);

    //--------------------------------//
    // Opzioni
    //
    string opz_dir = "../data/RR/";
    //string fname_opz = opz_dir + "opz-rr-EarthMars_test.yaml"; 
    string fname_opz = opz_dir + "opz-rr-PO.yaml"; 

    C_prb_RR_HMS prb(fname_opz);
    
    //--------------------------------//
    // inizializzazione
    //
    vector<int> nSeg_vector = prb.nSeg_vector;
    vector<double> tf_vector = prb.tfin_vector;
    int nLeg = prb.opz.nLeg;
    int nVars = prb.nVars;
    int opt_n = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg);
    //int opt_n = nLeg * (nSeg + 1) * prb.nVars;
    double *Xguess = new double[opt_n];
    for (int i = 0; i < opt_n; i++)
        Xguess[i] = 0;

    //--------------------------------//
    // Guess
    // 

    prb.init_lin(Xguess);
    string fname;
    if (prb.opz.firstguess_folder.compare("d") == 0)
        fname = "../results/HMS-temp/opt_X.dat";
    else
        fname = "../dbSOL/" + prb.opz.firstguess_folder + "/opt_X.dat";
    prb.load_sol(fname, Xguess);

    // for (int i = 0; i < opt_n; i++){
    //     for(int j=0; j<prb.nVars; j++){
    //         cout << Xguess[i*prb.nVars + j] << " ";
    //     }
    //     cout << endl;
    // }
    // cin.get();
    
    //--------------------------------//
    // Setup
    // 
    prb.setup_worph(Xguess);

    int MaxIter = prb.opz.MaxIter;

    time_t tstart, tend; 

    if (MaxIter>0){
        prb.par.MaxIter = MaxIter;
        tstart = time(0);
        prb.solve();
        tend = time(0); 
    }
    else{
        prb.opt.X = Xguess;
    }
    //--------------------------//
    //      Post process        //
    //--------------------------//
    string cmd;
    string output_dir = "../results/HMS-temp/";

    cmd = string("rm -r " + output_dir + "*"); system(cmd.c_str()); 
    system(string("mkdir -p " + output_dir).c_str());

    string fname_savings = "../results/HMS-temp/opt_X.dat";
    string out_opz("../results/HMS-temp/summary.yaml");
    double F = 0;
    double FLeg;

    int Nsub = prb.Nsub;

    MatrixXd P_kmf;
    VectorXd r_k, v_km;
    Vector3d r_k1, v_k1m, v_kp;
    MatrixXd P_km(6, 6), P_k1m(6, 6), P_kp(6, 6), Qd_k(6, 6), RR_k, cov_DV_i;
    VectorXd r_kt(2);

    VectorXd tspan;

    vector<double> v_DVnorm, v_DVstd, v_DVnormC, v_DVstdC;
    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "output";
        out << YAML::BeginMap;
    string namePm;  

    for(int iLeg = 0; iLeg < nLeg; iLeg++) {
        int nSeg = nSeg_vector[iLeg];
        cmd = "mkdir -p ../results/HMS-temp/sommario/"; system(cmd.c_str());
        cmd = "mkdir -p ../results/HMS-temp/Traj/"; system(cmd.c_str());
        cmd = "mkdir -p ../results/HMS-temp/DVs/"; system(cmd.c_str());
        cmd = "mkdir -p ../results/HMS-temp/MPm/"; system(cmd.c_str());
        cmd = "mkdir -p ../results/HMS-temp/MPp/"; system(cmd.c_str());
        cmd = "mkdir -p ../results/HMS-temp/PE/"; system(cmd.c_str());
        ofstream out_sommario("../results/HMS-temp/sommario/sommario" + to_string(iLeg + 1) + ".txt");
        string foldername = "../results/HMS-temp/Traj/Traj" + to_string(iLeg + 1) + ".dat";
        ofstream out_DVs("../results/HMS-temp/DVs/DVs" + to_string(iLeg + 1) + ".dat");
        ofstream out_MPm("../results/HMS-temp/MPm/MPm" + to_string(iLeg + 1) + ".dat");
        ofstream out_MPp("../results/HMS-temp/MPp/MPp" + to_string(iLeg + 1) + ".dat");
        ofstream out_covEellipsesXY("../results/HMS-temp/PE/PEXY" + to_string(iLeg + 1) + ".dat");
        ofstream out_covEellipsesXZ("../results/HMS-temp/PE/PEXZ" + to_string(iLeg + 1) + ".dat");
        ofstream out_covEellipsesYZ("../results/HMS-temp/PE/PEYZ" + to_string(iLeg + 1) + ".dat");

        ofstream out_sol(foldername);

        if (iLeg == 0){
            r_k = prb.r0;
            v_km = prb.v0;
            P_km = prb.P0;
        }
        else{
            r_k = prb.rRV_vector[iLeg - 1];
            v_km = v_kp;
            P_km = P_kp;
        }
        Qd_k = prb.Qd_k;

        tspan = VectorXd::LinSpaced(nSeg + 1, prb.t0, tf_vector[iLeg]);

        C_CR3BPMultiArc traj;
        // traj.add_CR3BPArc(C_CR3BPArc(r_k, v_km, prb.opz.mu, tspan(0), 0));
        
        for (int k = 0; k < nSeg + 1; k++){
            int i0 = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + k*prb.nVars;
            //int i0 = iLeg*(nSeg+1)*nVars + k*prb.nVars;
            Eigen::Map<Vector3d>  r_k(prb.opt.X + i0, 3);
            Eigen::Map<Vector3d> v_km(prb.opt.X + i0+3, 3);
            if ((prb.opz.v_RV_free)&&(k == 0)&&(iLeg > 0)){v_km = v_kp;}
            Eigen::Map<Vector3d> DV_i(prb.opt.X + i0+6, 3);
            Eigen::Map<MatrixXd> KK_i(prb.opt.X + i0+6+3, 3, 6);
            MatrixXd KKaug_i(6, 6); KKaug_i << MatrixXd::Zero(3, 6), KK_i;

            // Apply DV
            v_kp = v_km + DV_i;
            P_kp = (eye6 + KKaug_i) * P_km * (eye6 + KKaug_i).transpose();
            cov_DV_i = KK_i*P_km*KK_i.transpose();   

            double DVnorm_i = DV_i.norm(); 
            double DVstd_i;

            if (prb.opz.DVstd_model == 0) 
                DVstd_i = sqrt(cov_DV_i.trace());                     // sqrt(tr(DVstd))
            else{
                DVstd_i = sqrt(cov_DV_i.eigenvalues().real().maxCoeff());      // max autovalore
                // cout << cov_DV_i.eigenvalues() << "\t" << DVstd_i << endl;
            }

            v_DVnorm.push_back(DVnorm_i);
            v_DVstd.push_back(DVstd_i);
            v_DVnormC.push_back(DVnorm_i);
            v_DVstdC.push_back(DVstd_i);

            out_DVs << fixed << setw(12) << setprecision(6) << tspan[k];
            for(int ii=0; ii<3; ii++) out_DVs << scientific << setw(20) << setprecision(6) << r_k[ii];
            for(int ii=0; ii<3; ii++) out_DVs << scientific << setw(20) << setprecision(6) << DV_i[ii]; 
            out_DVs << endl;

            print_covariance_matrix(P_km, out_MPm);
            print_covariance_matrix(P_kp, out_MPp);
            out_MPm << endl;
            out_MPp << endl;

            // print_covariance(r_k, v_kp, P_km, out_covEellipses);
            print_covariance(r_k, remove_column_row(P_kp, 2, 2), out_covEellipsesXY);
            print_covariance(r_k(seq(1, 2)), remove_column_row(P_kp, 0, 0), out_covEellipsesYZ);
            r_kt << r_k(0), r_k(2); 
            print_covariance(r_kt, remove_column_row(P_kp, 1, 1), out_covEellipsesXZ);

            out_sommario << fixed << setw(2) << setprecision(8) << k << scientific << setw(2) << setprecision(6) << P_km.diagonal().transpose() << endl;
            out_sommario << fixed << setw(2) << setprecision(8) << k << scientific << setw(2) << setprecision(6) << P_kp.diagonal().transpose() << endl;
            out_sommario << endl;

            VectorXd Pm_vec = P_km.diagonal().transpose();
            std::vector<double> Pm(Pm_vec.data(), Pm_vec.data() + Pm_vec.size());

            namePm = "P_" + std::to_string(iLeg + 1) + ", " + std::to_string(k) + "-";
                out << YAML::Key << namePm << YAML::Value << Pm;
                out << YAML::Newline << YAML::Newline;

            if (k != nSeg){
                // propagate
                double dt = tspan(k + 1) - tspan(k);
                
                // Vector3d r_k1_hat, v_k1m_hat;  MatrixXd P_k1m_hat;
                propagate_UT(r_k, v_kp, P_kp, dt, prb.opz.mu, Qd_k, Nsub, r_k1, v_k1m, P_k1m);

                traj.add_CR3BPArc(C_CR3BPArc(r_k, v_kp, prb.opz.mu, tspan(k), dt));
            }
            else{
                r_k1 = r_k;
                v_k1m = v_kp;
                P_k1m = P_kp;
                P_kmf = P_km;

                // traj.add_CR3BPArc(C_CR3BPArc(r_k, v_kp, prb.opz.mu, tspan(k), 0));
            }

            // update
            P_km = P_k1m;
        }

        out_sommario << P_kmf.block(0,0,3,3).eigenvalues().real().transpose() << endl; // P_km instead of P_kmf?
        out_sommario << P_kmf.block(3,3,3,3).eigenvalues().real().transpose() << endl; // P_km instead of P_kmf?

        out_sommario << endl;

        for (int k=0; k<nSeg+1; k++){
            out_sommario << fixed << setw(2) << setprecision(8) << k 
                         << scientific << setw(18) << setprecision(10) << v_DVnorm[k] << "\t"
                         << scientific << setw(18) << setprecision(10) << v_DVstd[k] << "\t" 
                         << scientific << setw(18) << setprecision(10) << v_DVnorm[k] + prb.kstd*v_DVstd[k] << "\t" << endl;
        }
        FLeg = accumulate(v_DVnorm.begin(), v_DVnorm.end(), 0.) + prb.kstd*(accumulate(v_DVstd.begin(), v_DVstd.end(), 0.));
        F += FLeg;
        out_sommario << endl;  
        out_sommario << "FLeg = " << FLeg << endl;

        v_DVnorm.clear(); v_DVstd.clear();
        traj.print(out_sol, Nsub);
        if (iLeg == nLeg - 1){
            out_sommario << endl;  
            out_sommario << "F = " << F << endl;
            out_sommario << "F = " << fixed << F*1e3*prb.vconv << " [m/s]" << endl;
        }

        out_sommario.close();
        out_sol.close();
        out_DVs.close();
        out_MPm.close();
        out_MPp.close();
        out_covEellipsesXY.close();
        out_covEellipsesXZ.close();
        out_covEellipsesYZ.close();
    }
        out << YAML::Newline;
        out << YAML::Key << "DVnorm" << YAML::Value << v_DVnormC;
        out << YAML::Newline << YAML::Newline;

        out << YAML::Key << "DVstd" << YAML::Value << v_DVstdC;
        out << YAML::Newline << YAML::Newline;

        out << YAML::Key << "Fde" << YAML::Value << accumulate(v_DVnormC.begin(), v_DVnormC.end(), 0.) << YAML::Comment("Total deterministic DV");
        out << YAML::Key << "Fst" << YAML::Value << accumulate(v_DVstdC.begin(), v_DVstdC.end(), 0.) << YAML::Comment("Total stochastic DV");
        out << YAML::Key << "F_adi" << YAML::Value << F << YAML::Comment("Total DV");
        out << YAML::Key << "F" << YAML::Value << F*1e3*prb.vconv << YAML::Comment(" [m/s] Total DV");
    out << YAML::EndMap;

    // system("cd ..; gnuplot plot-files/plot_traj_RR_HMS.plt");

    // Save YP to Disk
    prb.save_sol(fname_savings);
    prb.opz.emit(out_opz);

    ofstream file_out(out_opz, ios::app);
    file_out << out.c_str() << endl;
    file_out.close();
    
    system("cd ..; gnuplot > /dev/null 2>&1 plot-files/plot_traj_RR_HMS.plt");

    string ans = "S";
    //cout << "Copy solution in dbSOL? (S/n)" << endl;
    //cin >> ans;
    if (ans.compare("S") == 0 || ans.compare("s") == 0){
            ofstream of("../results/HMS-temp/Traj.dat", std::ios_base::binary);

        for (int i = 1; i <= nLeg; i++) {
            string filename = "../results/HMS-temp/Traj/Traj" + to_string(i) + ".dat";
            ifstream infile(filename, ios::in | ios::binary);

            if (!infile) {
                cerr << "Error: Failed to open file " << filename << endl;
            }

            of << infile.rdbuf();

            infile.close();
        }

        of.close();

        const char* folder;
        string f = "/home/spectral01/MARMO/DESTINY+/PO/dbSOL/" + prb.opz.output_folder;
        folder = f.c_str();
        struct stat sb;

        if (stat(folder, &sb) == 0 && S_ISDIR(sb.st_mode)){
            cmd = string("rm -r ../dbSOL/" + prb.opz.output_folder); system(cmd.c_str());
        }
        cmd = "mkdir -p ../dbSOL/" + prb.opz.output_folder; system(cmd.c_str());
        cmd = string("cp -RT ") + "../results/HMS-temp" + " ../dbSOL/"+prb.opz.output_folder;  system(cmd.c_str()); 
        // cmd = "cp -RT " + output_dir + " ../dbSOL/"+prb.opz.output_folder;  system(cmd.c_str()); cout << cmd << endl;
    }

    cout << "It took " << difftime(tend, tstart) << " second(s) to solve the ROCP."<< endl;
}