#include "../include/C_prb_RR_HMS.h"
#include<string>
#include<sstream>
#include <ctime>

 void C_prb_RR_HMS::init_lin(double *Xguess){
    int krev = 0;

    // Default init: X=0
    for(int i = 0; i < nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); i++) 
        Xguess[i] = 0.;

    for(int iLeg = 0; iLeg < nLeg; iLeg++){
        int nSeg = nSeg_vector[iLeg];
        double tfin = tf_vector[iLeg];

        double r0_norm, th0, rf_norm, thf;
        Vector3d r0K, v0K, rfK, vfK;
    
        if(iLeg == 0){
            propagateKEP_U(r0_RV_vector[0], v0_RV_vector[0], tfin, 1, rfK, vfK);

            r0_norm = r0.norm();
            th0 = atan2(r0(1), r0(0));
            rf_norm = rfK.norm();
            thf = atan2(rfK(1), rfK(0)) + M_PI*2.*krev;
        }
        else{
            propagateKEP_U(r0_RV_vector[iLeg - 1], v0_RV_vector[iLeg - 1], accumulate(tf_vector.begin(), tf_vector.begin() + iLeg, 0.), 1, r0K, v0K);
            propagateKEP_U(r0_RV_vector[iLeg], v0_RV_vector[iLeg], accumulate(tf_vector.begin(), tf_vector.begin() + iLeg + 1, 0.), 1, rfK, vfK);

            r0_norm = r0K.norm();
            th0 = atan2(r0K(1), r0K(0));
            rf_norm = rfK.norm();
            thf = atan2(rfK(1), rfK(0)) + M_PI*2.*krev;
        }
        VectorXd v_r = VectorXd::LinSpaced(nSeg + 1, r0_norm, rf_norm);
        VectorXd v_th = VectorXd::LinSpaced(nSeg + 1, th0, thf);
        VectorXd v_vr = VectorXd::Zero(nSeg+1);
        VectorXd v_vt(nSeg + 1);
        for (int k = 0; k < nSeg + 1; k++)
            v_vt[k] = sqrt(1./v_r[k]);
        
        MatrixXd stateGuess(7, nSeg + 1);

        // r, th, vr, vt --> X,Y,VX,VY
        for (int k = 0; k < nSeg + 1; k++){
            double r = v_r[k];
            double th = v_th[k];
            double vr = v_vr[k];
            double vt = v_vt[k];

            double x = r*cos(th);
            double y = r*sin(th);
            double z = 0;
            double vx = vr*cos(th) - vt*sin(th);
            double vy = vr*sin(th) + vt*cos(th);
            double vz = 0;

            if(k < nSeg){
                stateGuess(0, k) = 0.5*tfin/nSeg;
            }
            else{
                stateGuess(0, k) = 0;
            }
            stateGuess(1, k) = x;
            stateGuess(2, k) = y;
            stateGuess(3, k) = z;
            stateGuess(4, k) = vx;
            stateGuess(5, k) = vy;
            stateGuess(6, k) = vz;     
        }

        int iX = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg);
        for (int k = 0; k < nSeg + 1; k++){
            int i0 = nVars*k;

            // state
            for (int j = 0; j < 6; j++)
                Xguess[iX + i0 + j] = stateGuess(j, k);
                        
            // DV [null]            

            // KK [null]
        }
    }
};

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

    cout << "rconv = " << rconv << endl;
    cout << "vconv = " << vconv << endl;
    cout << "tconv = " << tconv << endl;
    cout << "aconv = " << aconv << endl;

    tf_vector = opz.tfin_vector;
    nSeg_vector = opz.nSeg_vector;
    r0_RV_vector = opz.r0_RV_vector;
    v0_RV_vector = opz.v0_RV_vector;
    r0 = opz.r0;
    v0 = opz.v0;
    rf = opz.rf;
    vf = opz.vf;

    P0 = MatrixXd::Zero(6, 6);
    P0.diagonal() << opz.sigma2_r0, opz.sigma2_r0, opz.sigma2_r0, opz.sigma2_v0, opz.sigma2_v0, opz.sigma2_v0;
    
    Qd_k = MatrixXd::Zero(Pdim, Pdim); 
    if(opz.Qd_level == 0)
        Qd_k = MatrixXd::Zero(Pdim, Pdim); 
    else {
        X = VectorXd::Zero(Pdim);
        string Qd_level_array[5] = {"Qds", "Qdm", "Qdl", "Qdxl", "Qdxxl"};
        opt_Qd.open ("../src/" + Qd_level_array[opz.Qd_level - 1] + ".dat", ios::in); 
        for (int i = 0; i < 6; i++){
            getline(opt_Qd, line);
            stringstream ss(line); 
            for (int j = 0; j < 6; j++)
            {
                ss >> el;
                X(j) = stod(el);
            }
            Qd_k.row(i) = X;
        }
        opt_Qd.close();
        }

    const double beta = 0.05; // Pr{DV > DVmax} < 1 - beta
    kstd = sqrt(2*log(1./beta) + sqrt(3.)); // circa 4.18 @ beta=0.05

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
    if (status == DataError || status == InitError){
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
    opt.n = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); // {DT, rm, vm, DVx, DVy, DVz, K} x (N+1)
    
    opt.m = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.);     //continuità r,v punti interni + covarianza   
    opt.m += nRes*nLeg;              // vincoli all'arrivo 
    if(opz.E_Pf_constraint) opt.m += 6*nLeg;

    if(opz.E_DV_cstr) opt.m += accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg; //ciascun DV è vincolato

    if(opz.Fixed_ToF_Leg) opt.m += nLeg; // ToF of each leg is fixed

    opt.m += 2*opz.FB_Legs.size();        // FB constraints (norm(v_inf-) = norm(v_inf+), r_p >= r_p_min)

    if(opz.Limited_ToF) opt.m += 1; // Cumulative ToF

    // All derivatives for this problem have a sparse structure, so  set the amount of nonzeros here
    wsp.DF.nnz = WorhpMatrix_Init_Dense;  
    wsp.DG.nnz = WorhpMatrix_Init_Dense;  
    wsp.HM.nnz = WorhpMatrix_Init_Dense;  

    WorhpInit(&opt, &wsp, &par, &cnt);
    if (cnt.status != FirstCall){
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
    int nSeg, i0;
    double tfin, Np1L, LDV, UDV, LKK, UKK;
    bool FB;
    for(int iLeg = 0; iLeg < nLeg; iLeg++){
        nSeg = nSeg_vector[iLeg];
        tfin = tf_vector[iLeg];
        Np1L = accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg;

        FB = false;
        FB = find(opz.FB_Legs.begin(), opz.FB_Legs.end(), iLeg) != end(opz.FB_Legs);

        for(int iSeg = 0; iSeg < nSeg + 1; iSeg++){
            i0 = nVars*Np1L + iSeg*nVars;

            if(iSeg < nSeg){ //DT
                if(opz.Equal_Segment_ToF){
                    opt.XL[i0] = tfin/nSeg;
                    opt.XU[i0] = tfin/nSeg;}
                else{
                    opt.XL[i0] = 0.25*tfin/nSeg;
                    opt.XU[i0] = tfin;}
            }
            else{
                opt.XL[i0] = .0;
                opt.XU[i0] = .0;
            }

            for(int j = 0; j < 3; j++){ //pos
                opt.XL[i0 + 1 + j] = -4.0;
                opt.XU[i0 + 1 + j] = +4.0;
            }
            for(int j = 0; j < 3; j++){ //vel
                opt.XL[i0 + 1 + 3 + j] = -2.0;
                opt.XU[i0 + 1 + 3 + j] = +2.0;
            }

            if((iLeg != 0)&&(iSeg == 0)){
                for (int j = 0; j < 3; j++){
                opt.XL[i0 + 1 + j] = opt.XU[i0 + 1 + j] = 0;
                opt.XL[i0 + 1 + 3 + j] = opt.XU[i0 + 1 + 3 + j] = 0;
                }  
            }
            
            // DV
            if((iSeg == nSeg)||((iSeg == 0)&&(!FB)&&(iLeg != 0))){
                    LDV = UDV = .0;
            }
            else if((iSeg == 0)&&(FB||(iLeg == 0))){
                LDV = -2;
                UDV = +2;
            }
            else{
                LDV = -opz.DVtot_single_max;
                UDV = +opz.DVtot_single_max;
            }
            for(int j = 0; j < 3; j++){ 
                if(opz.E_DV_cstr){   
                    opt.XL[i0 + 1 + 6 + j] = LDV;
                    opt.XU[i0 + 1 + 6 + j] = UDV;
                }
                else{
                    opt.XL[i0 + 1 + 6 + j] = -1.0;
                    opt.XU[i0 + 1 + 6 + j] = +1.0;
                }
            }

            // KK
            if(opz.E_Pf_constraint){
                if((iSeg == nSeg)||(iSeg == 0)){
                    LKK = UKK = .0;
                }
                else{
                    LKK = -200.0;
                    UKK = +200.0;
                }
                //if(iLeg > 0) LKK = UKK = 0;
                for(int iRow = 0; iRow < 3; iRow++){       
                    for(int iCol = 0; iCol < 6; iCol++){
                        opt.XL[i0 + 1 + 6 + 3 + iRow*6 + iCol] = LKK;
                        opt.XU[i0 + 1 + 6 + 3 + iRow*6 + iCol] = UKK;
                    }
                }
            }
            else{
                for(int iRow = 0; iRow < 3; iRow++){       
                    for(int iCol = 0; iCol < 6; iCol++){
                        opt.XL[i0 + 1 + 6 + 3 + iRow*6 + iCol] = opt.XU[i0 + 6 + 3 + 1 + iRow*6 + iCol] = 0;
                    }
                }
            }
        }

        // Condizioni iniziali assegnate
        if(iLeg == 0){                          
            for(int j = 0; j < 3; j++){
                opt.XL[1 + j] = opt.XU[1 + j] = r0(j);
                opt.XL[1 + 3 + j] = opt.XU[1 + 3 + j] = v0(j);
            }  
        }

        // Vincoli: Residui + Condizioni al contorno
        // residuo r & v
        for(int iSeg = 0; iSeg < nSeg; iSeg++){
            int iG0 = (nRes + 1)*Np1L + iSeg*(nRes); // r, v, sigma2_rf, sigma2_vf, DV
            if(opz.E_Pf_constraint) iG0 += 6*iLeg;
            if(opz.Fixed_ToF_Leg) iG0 += iLeg;

            //residuo_rv
            for(int j = 0; j < 6; j++){  
                opt.GL[iG0 + j] = opt.GU[iG0 + j] = 0; 
            }
        }

        // condizioni contorno finale
        int iG0 = (nRes + 1)*Np1L + nSeg*(nRes);
        if(opz.E_Pf_constraint) iG0 += 6*iLeg;
        if(opz.Fixed_ToF_Leg) iG0 += iLeg;

        for(int i = 0; i < 6; i++){
            opt.GL[iG0 + i] = opt.GU[iG0 + i] = 0; 
        }

        // covarianza di posizione finale: Pr - Prdes < 0
        if(opz.E_Pf_constraint){
            for(int j = 0; j < 6; j++){
                opt.GL[iG0 + 6 + j] = -par.Infty;   
                opt.GU[iG0 + 6 + j] = 0;
            }
        } 
        // DVcstr: DV[k] - DV_max < 0
        if (opz.E_DV_cstr){
            int nPf_cstr = 0;
            if(opz.E_Pf_constraint) nPf_cstr = 6;
            for(int j = 0; j < nSeg + 1; j++){
                opt.GL[iG0 + 6 + nPf_cstr + j] = -par.Infty; 
                opt.GU[iG0 + 6 + nPf_cstr + j] = 0;
            }
        } 

        // DT constraints
        int nPf_cstr = 0;
        if(opz.E_Pf_constraint) nPf_cstr = 6;
        if(opz.Fixed_ToF_Leg){
            opt.GL[iG0 + 6 + nPf_cstr + nSeg + 1] = opt.GU[iG0 + 6 + nPf_cstr + nSeg + 1] = 0;
        }
    }
    int iG0 = (nRes + 1)*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); 
    if(opz.E_Pf_constraint) iG0 += 6*nLeg;
    if(opz.Fixed_ToF_Leg) iG0 += nLeg;
    for(int i; i < opz.FB_Legs.size(); i++){
        opt.GL[iG0 + 2*i] = opt.GU[iG0 + 2*i] = 0; // norm(v_inf-) = norm(v_inf+)
        opt.GL[iG0 + 2*i + 1] = 0; // r_p >= r_p_min
        opt.GU[iG0 + 2*i + 1] = opz.SoI_R; // SoI radius
    }

    if(opz.Limited_ToF){ // ToF - Max_ToF <= 0
        opt.GL[opt.m - 1] = -par.Infty;
        opt.GU[opt.m - 1] = 0;
    }
    par.FGtogether = true;
    par.UserDF = false;
    par.UserDG = false;
    par.UserHM = false;

    // for(int i = 0; i < opt.n; i++){
    //         cout << "XL[" << i << "] = " << opt.XL[i] << endl;
    // }
    // for(int i = 0; i < opt.n; i++){
    //     cout << "XU[" << i << "] = " << opt.XU[i] << endl;
    // }
    cin.get();
    return 0;
}

void C_prb_RR_HMS::custom_iteration_process(double *X_current){
    return;
}

void C_prb_RR_HMS::evalFG(double *X, double &F, double *G, double ScaleObj){
    vector<Vector3d> v_infm_vec, v_infp_vec;
    vector<double> r_p_vec;
    MatrixXd eye6x6 = MatrixXd::Identity(6, 6);
    C_UT ut(6);
    MatrixXd P_k1m(6, 6), P_km(6, 6), P_kp(6, 6), covDV(3, 3), KK_k(3, 6), KKaug_k(6, 6), P_k1m_hat(6, 6), KK_k1(3, 6), P_FB(6, 6);
    Vector3d r_km, v_km, r_k1m, v_k1m, v_kp, v_k1p, DV_k, DV_k1, r_k1_hat, v_k1m_hat, Pf_ev_r, Pf_ev_v, v_infp, v_infm;

    F = 0;
    
    int nSeg, FB_c;
    double tfin, Np1L, dt, T, ToF, theta, e, a, v_inf_sqrd, r_p;
    ToF = 0;
    FB_c = 0;
    bool FB, FBm1;
    for(int iLeg = 0; iLeg < nLeg; iLeg++){
        nSeg = nSeg_vector[iLeg];
        tfin = tf_vector[iLeg];
        T = 0;
        vector<double> v_DVnorm(nSeg + 1); 
        vector<double> v_DVstd(nSeg + 1);

        Np1L = accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg;
        int iX = nVars*Np1L;

        dt = X[iX];
        r_km = Eigen::Map<Vector3d>(X + iX + 1, 3);  // Initial state is assigned by setting upper and lower boundaires in X equal to \tilde{x0} 
        v_km = Eigen::Map<Vector3d>(X + iX + 1 + 3, 3);    

        if(iLeg == 0){
            P_km = P0;  
        }        
        else{
            P_km = P_kp;
            if(opz.v_RV_free){
                r_km = r_k1m;
                v_km = v_kp;
            }
        }                          
        DV_k = Eigen::Map<Vector3d>(X + iX + 1 + 6, 3);      
        KK_k = Eigen::Map<MatrixXd>(X + iX + 1 + 6 + 3, 3, 6); 
        KKaug_k << MatrixXd::Zero(3, 6), KK_k; 

        // v_infp (FB already "performed" in previous cycle)
        FB = false;
        FB = find(opz.FB_Legs.begin(), opz.FB_Legs.end(), iLeg) != end(opz.FB_Legs);
        if(FB){
            FB_c++;
            propagateKEP_U(r0_RV_vector[iLeg - 1], v0_RV_vector[iLeg - 1], ToF, 1, r_G, v_G);
            v_infp = DV_k - v_G + v_kp; // DV - v_p + v_1
            v_infp_vec.push_back(v_infp);
            // Propagate_P_FB(P_km, v_infp_vec[FB_c - 1], v_infp, 929000, opz.mu_FB, P_FB);
            // P_km = P_FB;
        }

        // apply DV
        v_kp = v_km + DV_k;
        P_kp = (eye6x6 + KKaug_k) * P_km * (eye6x6 + KKaug_k).transpose();
        covDV = KK_k*P_km*KK_k.transpose();

        if(FB||(iLeg == 0)){
            v_DVnorm[0] = 0; 
            v_DVstd[0] = 0;
        }
        else{
            v_DVnorm[0] = DV_k.norm(); 
            v_DVstd[0] = sqrt(covDV.eigenvalues().real().maxCoeff()); //max autovalore
        }

        for(int iSeg = 0; iSeg < nSeg; iSeg++){
            T += dt;
            // ottenuti per propagazione
            propagate_kepler_UT(r_km, v_kp, P_kp, dt, 1., Qd_k, r_k1_hat, v_k1m_hat, P_k1m_hat);

            int i0 = iX + (iSeg + 1)*nVars;
            dt = X[i0];
            r_k1m = Eigen::Map<Vector3d>(X + i0 + 1, 3);  
            v_k1m = Eigen::Map<Vector3d>(X + i0 + 1 + 3, 3);                                       
            P_k1m = P_k1m_hat;                                  // Hybrid: Covarianza propagata
            DV_k1 = Eigen::Map<Vector3d>(X + i0 + 1 + 6, 3);        
            KK_k1 = Eigen::Map<MatrixXd>(X + i0 + 1 + 6 + 3, 3, 6);        
            
            int iG0 = (nRes + 1)*Np1L + iSeg*(nRes);
            if(opz.E_Pf_constraint) iG0 += 6*iLeg;
            if(opz.Fixed_ToF_Leg) iG0 += iLeg;

            for(int j = 0; j < 3; j++) G[iG0 + j] = r_k1m[j] - r_k1_hat[j];
            for(int j = 0; j < 3; j++) G[iG0 + 3 + j] = v_k1m[j] - v_k1m_hat[j];

            // apply DV
            v_k1p = v_k1m + DV_k1;
            MatrixXd KKaug_k1(6, 6); KKaug_k1 << MatrixXd::Zero(3, 6), KK_k1;  
            MatrixXd P_k1p = (eye6x6 + KKaug_k1) * P_k1m_hat * (eye6x6 + KKaug_k1).transpose();  
            MatrixXd cov_DV_k1 = KK_k1*P_k1m*KK_k1.transpose();          
            

            v_DVnorm[iSeg + 1] = DV_k1.norm(); 
            v_DVstd[iSeg + 1] = sqrt(cov_DV_k1.eigenvalues().real().maxCoeff());      // max autovalore

            // update
            r_km = r_k1m; // note that: r_km = r_kp
            v_kp = v_k1p;
            P_kp = P_k1p;
        }
        ToF += T;

        // v_infm (FB to be "performed" now)
        FBm1 = false;
        FBm1 = find(opz.FB_Legs.begin(), opz.FB_Legs.end(), iLeg + 1) != end(opz.FB_Legs);
        if(FBm1){
            propagateKEP_U(r0_RV_vector[iLeg], v0_RV_vector[iLeg], ToF, 1, r_G, v_G);
            v_infm = v_kp - v_G; // v_1 - v_p
            v_infm_vec.push_back(v_infm);
        }

        // obj function
        double DVtot = accumulate(v_DVnorm.begin(), v_DVnorm.end(), 0.);
        double DVstd_tot = accumulate(v_DVstd.begin(), v_DVstd.end(), 0.);
        //     F = DVtot;
        F += DVtot + DVstd_tot*kstd;   

        propagateKEP_U(r0_RV_vector[iLeg], v0_RV_vector[iLeg], ToF, 1, r_G, v_G);

        // 1.1) Stato medio al tempo finale
        int iG0 = (nRes + 1)*Np1L + nSeg*(nRes);
        if(opz.E_Pf_constraint) iG0 += 6*iLeg;
        if(opz.Fixed_ToF_Leg) iG0 += iLeg;

        if(iLeg == nLeg - 1){
            for(int j = 0; j < 3; j++) G[iG0 + j] = r_km[j] - r_G[j];
            for(int j = 0; j < 3; j++) G[iG0 + 3 + j] = v_kp[j] - v_G[j];
        }
        else{
            for(int j = 0; j < 3; j++) G[iG0 + j] = r_km[j] - r_G[j];
            for(int j = 0; j < 3; j++) G[iG0 + 3 + j] = 0;
        }

        // 1.2) Covarianza al tempo finale
        if(opz.E_Pf_constraint){
            Pf_ev_r = P_kp.block(0,0,3,3).eigenvalues().real();
            Pf_ev_v = P_kp.block(3,3,3,3).eigenvalues().real();
            
            if(iLeg == nLeg - 1){
                sigma2_r_G = opz.sigma2_rf;
                sigma2_v_G = opz.sigma2_vf;
                }
            else{
                sigma2_r_G = opz.sigma2_rRV;
                sigma2_v_G = opz.sigma2_vRV;
            }
             
            for(int j = 0; j < 3; j++){
                //G[iG0 + 6 + j] = P_kp(j,j) - sigma2_r_G; 
                G[iG0 + 6 + j] = opz.sf*(Pf_ev_r(j) - sigma2_r_G);
            }
            for(int j = 3; j < 6; j++){
                //G[iG0 + 6 + j] = P_kp(j,j) - sigma2_v_G; 
                G[iG0 + 6 + j] = opz.sf*(Pf_ev_v(j-3) - sigma2_v_G); 
            }
        }

        // 1.3) DV_cstr
        if(opz.E_DV_cstr){
            int nPf_cstr = 0;
            if(opz.E_Pf_constraint) nPf_cstr = 6;{
                for(int k = 0; k < nSeg + 1; k++){
                    G[iG0 + 6 + nPf_cstr + k] = v_DVnorm[k] + kstd*v_DVstd[k] - opz.DVtot_single_max; 
                }
            }
        }

        int nPf_cstr = 0;
        if(opz.E_Pf_constraint) nPf_cstr = 6;
        if(opz.Fixed_ToF_Leg){
            G[iG0 + 6 + nPf_cstr + nSeg + 1] = T - tfin;
        }
    }

    int iG0 = (nRes + 1)*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg);
    if(opz.E_Pf_constraint) iG0 += 6*nLeg;
    if(opz.Fixed_ToF_Leg) iG0 += nLeg;
    for(int i; i < opz.FB_Legs.size(); i++){
        v_inf_sqrd = v_infm_vec[i].squaredNorm();
        theta = acos((v_infm_vec[i]).dot(v_infp_vec[i])/v_inf_sqrd);
        r_p = opz.mu_FB*(1 - sin(theta/2))/(v_inf_sqrd*vconv*vconv*sin(theta/2));
        G[iG0 + 2*i] = v_infm_vec[i].norm() - v_infp_vec[i].norm(); // v_infm_vec[i].norm() - v_infp_vec[i].norm()
        G[iG0 + 2*i + 1] = r_p - opz.r_min; // r_p - opz.r_min
    }
    if(opz.Limited_ToF) G[opt.m - 1] = ToF - opz.Max_ToF;
    F *= ScaleObj;
}
 
void test_RR_HMS(){
    MatrixXd eye6 = MatrixXd::Identity(6, 6);

    // Opzioni
    string opz_dir = "../data/RR/";
    //string fname_opz = opz_dir + "opz-rr-EarthMars_test.yaml"; 
    string fname_opz = opz_dir + "opz-rr-AFC1.yaml"; 

    C_prb_RR_HMS prb(fname_opz);
    
    // inizializzazione
    vector<int> nSeg_vector = prb.nSeg_vector;
    vector<double> tf_vector = prb.tf_vector;
    int nLeg = prb.opz.nLeg;
    int nVars = prb.nVars;
    int opt_n = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg);
    //int opt_n = nLeg * (nSeg + 1) * prb.nVars;
    double *Xguess = new double[opt_n];
    for (int i = 0; i < opt_n; i++)
        Xguess[i] = 0;

    // Guess
    prb.init_lin(Xguess);
    string fname;
    if(prb.opz.firstguess_folder.compare("d") == 0)
        fname = "../results/HMS-temp/opt_X.dat";
    else
        fname = "../dbSOL/" + prb.opz.firstguess_folder + "/opt_X.dat";
    prb.load_sol(fname, Xguess);
    
    // Setup
    prb.setup_worph(Xguess);

    int MaxIter = prb.opz.Max_iter;

    time_t tstart, tend; 

    if(MaxIter > 0){
        prb.par.MaxIter = MaxIter;
        tstart = time(0);
        prb.solve();
        tend = time(0); 
    }
    else
        prb.opt.X = Xguess;

    //--------------------------//
    //      Post process        //
    //--------------------------//
    string cmd;
    string output_dir = "../results/HMS-temp/";

    cmd = string("rm -r " + output_dir + "*"); system(cmd.c_str()); 
    system(string("mkdir -p " + output_dir).c_str());

    ofstream trajEfile("../results/HMS-temp/trajEfile.dat");
    ofstream trajRV1file("../results/HMS-temp/trajRV1file.dat");
    ofstream trajRV3file("../results/HMS-temp/trajRV3file.dat");
    string fname_savings = "../results/HMS-temp/opt_X.dat";
    string out_opz("../results/HMS-temp/opz.yaml");
    double F, FLeg, dt, ToF, ToF_Leg, tfin, DVnorm_i, DVstd_i;
    int nSeg;
    F = 0;
    ToF = 0;

    Vector3d r_k, v_km, r_k1, v_k1m, v_kp, DV_i, v_infp, v_infm, r_G, v_G;
    MatrixXd P_km(6, 6), P_k1m(6, 6), P_kp(6, 6), Qd_k(6, 6), cov_DV_i(6, 6), KK_i(3, 6), KKaug_i(6, 6);
    vector<Vector3d> v_infm_vec, v_infp_vec;

    vector<double> v_DVnorm, v_DVstd, ToF_par;
    C_KeplerMultiArc trajE, trajRV1, trajRV3;
    bool FB, FBm1;

    for(int iLeg = 0; iLeg < nLeg; iLeg++){
        ToF_Leg = 0;
        nSeg = nSeg_vector[iLeg];
        tfin = tf_vector[iLeg];
        ofstream out_sommario(output_dir + "sommario" + to_string(iLeg + 1) + ".txt");
        string foldername = "../results/HMS-temp/fullsol" + to_string(iLeg + 1) + ".dat";
        ofstream out_DVs("../results/HMS-temp/DVs" + to_string(iLeg + 1) + ".dat");
        ofstream out_covEellipses("../results/HMS-temp/covEllipses" + to_string(iLeg + 1) + ".dat");
        ofstream out_sol(foldername);

        if(iLeg == 0){
            r_k = prb.r0;
            v_km = prb.v0;
            P_km = prb.P0;
        }
        else{
            r_k = r_k1;
            v_km = v_kp;
            P_km = P_kp;
        }
        Qd_k = prb.Qd_k;

        C_KeplerMultiArc traj;
        traj.add_keplerArc(C_KeplerArc(r_k, v_km, 0, 0));

        for(int iSeg = 0; iSeg < nSeg + 1; iSeg++){
            int i0 = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + iSeg*prb.nVars;
            dt = prb.opt.X[i0];
            r_k = Eigen::Map<Vector3d>(prb.opt.X + i0 + 1, 3);
            v_km = Eigen::Map<Vector3d>(prb.opt.X + i0 + 1 + 3, 3);
            if((prb.opz.v_RV_free)&&(iSeg == 0)&&(iLeg > 0)){v_km = v_kp; r_k = r_k1;}
            DV_i = Eigen::Map<Vector3d>(prb.opt.X + i0 + 1 + 6, 3);
            KK_i = Eigen::Map<MatrixXd>(prb.opt.X + i0 + 1 + 6 + 3, 3, 6);
            KKaug_i << MatrixXd::Zero(3, 6), KK_i;   

            // v_infp
            if(iSeg == 0){
                FB = false;
                FB = find(prb.opz.FB_Legs.begin(), prb.opz.FB_Legs.end(), iLeg) != end(prb.opz.FB_Legs);
                if(FB){
                    propagateKEP_U(prb.r0_RV_vector[iLeg - 1], prb.v0_RV_vector[iLeg - 1], ToF, 1, r_G, v_G);
                    v_infp = DV_i - v_G + v_kp;
                    v_infp_vec.push_back(v_infp);
                }
            }

            // Apply DV
            v_kp = v_km + DV_i;
            P_kp = (eye6 + KKaug_i) * P_km * (eye6 + KKaug_i).transpose();
            cov_DV_i = KK_i*P_km*KK_i.transpose();

            if((FB&&(iSeg == 0))||(iLeg == 0)&&(iSeg == 0)){
                DVnorm_i = 0; 
            }
            else{
                DVnorm_i = DV_i.norm(); 
            }
            double DVnorm_i = DV_i.norm(); 
            DVstd_i = sqrt(cov_DV_i.eigenvalues().real().maxCoeff());      // max autovalore
        
            v_DVnorm.push_back(DVnorm_i);
            v_DVstd.push_back(DVstd_i);

            out_DVs << fixed << setw(12) << setprecision(6) << ToF;
            for(int ii=0; ii<3; ii++) out_DVs << fixed <<setw(12) << setprecision(6) << r_k[ii];
            for(int ii=0; ii<3; ii++) out_DVs << fixed <<setw(12) << setprecision(6) << DV_i[ii]; 
            out_DVs << endl;

            print_covariance(r_k, v_kp, P_km, out_covEellipses);
            print_covariance(r_k, v_kp, P_kp, out_covEellipses);

            out_sommario << fixed << setw(2) << setprecision(8) << iSeg << scientific << setw(2) << setprecision(5) << P_km.diagonal().transpose() << endl;
            out_sommario << fixed << setw(2) << setprecision(8) << iSeg << scientific << setw(2) << setprecision(5) << P_kp.diagonal().transpose() << endl;
            out_sommario << endl;

            if (iSeg != nSeg){
                // propagate
                
                // Vector3d r_k1_hat, v_k1m_hat;  MatrixXd P_k1m_hat;
                propagate_kepler_UT(r_k, v_kp, P_kp, dt, 1., Qd_k, r_k1, v_k1m, P_k1m);
                traj.add_keplerArc(C_KeplerArc(r_k, v_kp, ToF, dt));
            }
            else{
                r_k1 = r_k;
                v_k1m = v_kp;
                P_k1m = P_kp;
                traj.add_keplerArc(C_KeplerArc(r_k, v_kp, ToF, 0));
            }
            // v_infm
            if(iSeg == nSeg){
                FBm1 = false;
                FBm1 = find(prb.opz.FB_Legs.begin(), prb.opz.FB_Legs.end(), iLeg + 1) != end(prb.opz.FB_Legs);
                if(FBm1){
                    propagateKEP_U(prb.r0_RV_vector[iLeg], prb.v0_RV_vector[iLeg], ToF, 1, r_G, v_G);
                    v_infm = v_kp - v_G;
                    v_infm_vec.push_back(v_infm);
                }
            }

            // update
            P_km = P_k1m;
            ToF += dt;
            ToF_Leg += dt;
        }
        ToF_par.push_back(ToF_Leg);

        out_sommario << P_km.block(0, 0, 3, 3).eigenvalues().real().transpose() << endl;
        out_sommario << P_km.block(3, 3, 3, 3).eigenvalues().real().transpose() << endl;

        out_sommario << endl;
        for (int k = 0; k < nSeg + 1; k++){
            out_sommario << fixed << setw(2) << setprecision(8) << k 
                << fixed << setw(18) << setprecision(10) << v_DVnorm[k] << "\t"
                << fixed << setw(18) << setprecision(10) << v_DVstd[k] << "\t" 
                << fixed << setw(18) << setprecision(10) << v_DVnorm[k] + prb.kstd*v_DVstd[k] << "\t" << endl;
        }
        FLeg = accumulate(v_DVnorm.begin(), v_DVnorm.end(), 0.) + prb.kstd*(accumulate(v_DVstd.begin(), v_DVstd.end(), 0.));
        F += FLeg;
        out_sommario << endl;  
        out_sommario << "FLeg = " << FLeg << endl;
        out_sommario << "ToF_Leg (actual) = " << ToF_Leg << endl;
        out_sommario << "ToF_Leg (design) = " << tfin << endl;

        v_DVnorm.clear(); v_DVstd.clear();
        traj.print(out_sol, 100);
        if (iLeg == nLeg - 1){
            out_sommario << endl;  
            out_sommario << "F = " << F << endl;

            out_sommario << "ToF (actual) = " << ToF << endl;
            out_sommario << "ToF (max) \t = " << prb.opz.Max_ToF << endl << endl;

            double v_inf_sqrd, theta, r_p;
            for(int i; i < prb.opz.FB_Legs.size(); i++){
                v_inf_sqrd = v_infm_vec[i].squaredNorm();
                theta = acos((v_infm_vec[i]).dot(v_infp_vec[i])/v_inf_sqrd);
                r_p = prb.opz.mu_FB*(1 - sin(theta/2))/(v_inf_sqrd*prb.opz.vconv*prb.opz.vconv*sin(theta/2));
                Propagate_P_FB(P_km, v_infm_vec[i], v_infp_vec[i], prb.opz.SoI_R/prb.opz.rconv, prb.opz.mu_FB/prb.opz.mu_primary, P_km);
                out_sommario << "FB #" << i + 1 << ":" << endl;
                out_sommario << "r_p = " << r_p << endl;
                out_sommario << "v_infm = " << v_infm_vec[i][0] << " " << v_infm_vec[i][1] << " " << v_infm_vec[i][2] << endl;
                out_sommario << "v_infp = " << v_infp_vec[i][0] << " " << v_infp_vec[i][1] << " " << v_infp_vec[i][2] << endl;
                out_sommario << "|v_infm| = " << v_infm_vec[i].norm() << endl;
                out_sommario << "|v_infp| = " << v_infp_vec[i].norm() << endl << endl;
            }
        }

        out_sommario.close();
        out_sol.close();
        out_DVs.close();
        out_covEellipses.close();
    }
    trajE.add_keplerArc(C_KeplerArc(prb.r0, prb.v0, 0, ToF));
    trajE.print(trajEfile, 200);
    trajEfile.close();

    trajRV1.add_keplerArc(C_KeplerArc(prb.r0_RV_vector[0], prb.v0_RV_vector[0], 0, accumulate(ToF_par.begin(), ToF_par.begin() + 1, 0.)));
    trajRV1.print(trajRV1file, 200);
    trajRV1file.close();

    trajRV3.add_keplerArc(C_KeplerArc(prb.r0_RV_vector[2], prb.v0_RV_vector[2], 0, accumulate(ToF_par.begin(), ToF_par.begin() + 3, 0.)));
    trajRV3.print(trajRV3file, 200);
    trajRV3file.close();

    // Save YP to Disk
    prb.save_sol(fname_savings);
    prb.opz.emit(out_opz);
    
    system("cd ..; gnuplot > /dev/null 2>&1 plot-files/plot_traj_RR_HMS.plt");

    string ans = "S";
    //cout << "Copy solution in dbSOL? (S/n)" << endl;
    //cin >> ans;
    string f;
    struct stat sb;

    f = "../dbSOL/" + prb.opz.output_folder;
    const char* folder = f.c_str();
    if(ans.compare("S") == 0 || ans.compare("s") == 0){
        ofstream of("../results/HMS-temp/fullsol.dat", ios_base::binary);
        for(int i; i < nLeg; i++){
            ifstream if_a("../results/HMS-temp/fullsol" + to_string(i + 1) + ".dat", ios_base::binary);
            of << if_a.rdbuf();
        }

        string cmd;
        if(stat(folder, &sb) == 0 && S_ISDIR(sb.st_mode)){
            cmd = string("cd ..; rm -r dbSOL/") + prb.opz.output_folder; system(cmd.c_str()); 
        }
        cmd = "mkdir -p ../dbSOL/" + prb.opz.output_folder; system(cmd.c_str()); //cout << cmd << endl;
        cmd = string("cp -RT ") + "../results/HMS-temp" + " ../dbSOL/" + prb.opz.output_folder;  system(cmd.c_str()); //cout << cmd << endl;
    }
    cout << "It took " << difftime(tend, tstart) << " second(s) to solve the ROCP."<< endl;
}