#include "../include/C_prb_RR_HMS.h"
#include<string>
#include<sstream>
#include <ctime>

 void C_prb_RR_HMS::init_lin(double *Xguess){
    int krev = 0;

    // Default init: X=0
    for (int i = 0; i < nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); i++) 
        Xguess[i] = 0.;

    for (int iLeg=0; iLeg < nLeg; iLeg++){
        int nSeg = nSeg_vector[iLeg];
        double tfin = tf_vector[iLeg];

        double r0_norm;
        double th0;

        double rf_norm;
        double thf;
        
        if (iLeg == 0){
            r0_norm = r0.norm();
            th0 = atan2(r0(1),r0(0));
            rf_norm = rRV_vector[0].norm();
            thf = atan2(rRV_vector[0](1), rRV_vector[0](0)) + M_PI*2.*krev;
        }
        else if(iLeg == nLeg - 1){
            r0_norm = rRV_vector[nLeg - 2].norm();
            th0 = atan2(rRV_vector[nLeg - 2](1), rRV_vector[nLeg - 2](0));
            rf_norm = rf_des.norm();
            thf = atan2(rf_des(1), rf_des(0)) + M_PI*2.*krev;
        }
        else{
            r0_norm = rRV_vector[iLeg - 1].norm();
            th0 = atan2(vRV_vector[iLeg - 1](1), vRV_vector[iLeg - 1](0));
            rf_norm = rRV_vector[iLeg].norm();
            thf = atan2(rRV_vector[iLeg](1), rRV_vector[iLeg](0)) + M_PI*2.*krev;
        }
        VectorXd v_r = VectorXd::LinSpaced(nSeg+1, r0_norm, rf_norm);
        VectorXd v_th = VectorXd::LinSpaced(nSeg+1, th0, thf);
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
    amrif = 1000; //kgs
    const double g0_dim = 9.81e-3; // km/s^2
    const double Isp_dim = 2000;    //s
    const double THR_dim = 0.5e-3; // kN
    THR = THR_dim / (amrif*aconv);  // nondimensional
    C = g0_dim * Isp_dim / vconv;
    
    //nSeg = opz.nSeg;
    nSeg1 = opz.nSeg1;
    nSeg2 = opz.nSeg2;
    nLeg = opz.nLeg;
    sf = opz.sf;

    cout << "rconv = " << rconv << endl;
    cout << "vconv = " << vconv << endl;
    cout << "tconv = " << tconv << endl;
    cout << "aconv = " << aconv << endl;

    t0 = 0;
    // tf1 = opz.tfin1;
    // tf2 = opz.tfin2;
    tf_vector = opz.tfin_vector;
    nSeg_vector = opz.nSeg_vector;
    rRV_vector = opz.rRV_vector;
    r0 = opz.r0;
    v0 = opz.v0;
    rf_des = opz.rf;
    vf_des = opz.vf;
    rRV = opz.rRV;
    vRV = opz.vRV;

    P0 = MatrixXd::Zero(6, 6);
    if (opz.planar)
        P0.diagonal() << opz.sigma_r0, opz.sigma_r0, 0., opz.sigma_v0, opz.sigma_v0, 0;
    else
        P0.diagonal() << opz.sigma_r0, opz.sigma_r0, opz.sigma_r0, opz.sigma_v0, opz.sigma_v0, opz.sigma_v0;
    
    Qd_k = MatrixXd::Zero(Pdim, Pdim); 
    if (opz.planar)
        Qd_k.diagonal() << opz.Qd_r, opz.Qd_r, 0., opz.Qd_v, opz.Qd_v, 0;
    else {
        //Qd_k.diagonal() << opz.Qd_r, opz.Qd_r, opz.Qd_r, opz.Qd_v, opz.Qd_v, opz.Qd_v;
        //Qd_k.diagonal() << opz.Qd_r, opz.Qd_r, 0., opz.Qd_v, opz.Qd_v, 0;
        if (opz.Qd_level == 0)
            Qd_k = MatrixXd::Zero(Pdim, Pdim); 
        else {
            X = VectorXd::Zero(Pdim);
            string Qd_level_array[5] = {"Qds", "Qdm", "Qdl", "Qdxl", "Qdxxl"};
            opt_Qd.open ("../src/" + Qd_level_array[opz.Qd_level - 1] + ".dat", ios::in); 
            for (int i = 0; i < 6; i++)
            {
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
    }
    
    RR_k = MatrixXd::Zero(Pdim, Pdim); 
    if (opz.planar)
        RR_k.diagonal() << opz.nav_sigma_r, opz.nav_sigma_r, 0., opz.nav_sigma_v, opz.nav_sigma_v, 0;
    else
        RR_k.diagonal() << opz.nav_sigma_r, opz.nav_sigma_r, opz.nav_sigma_r, opz.nav_sigma_v, opz.nav_sigma_v, opz.nav_sigma_v;

    // Probability ADIM
    if (opz.cov_propagation_mode == 0)
        prconv = 1.;
    else
        prconv = 1;

    P0 /= prconv;   
    Qd_k /= prconv;
    RR_k /= prconv;

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
    opt.n = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg); // {DT, rm, vm, DVx, DVy, DVz, K} x (N+1)
    
    opt.m = 6*nLeg;              // vincoli all'arrivo
    opt.m += nRes*accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.);     //continuità r,v punti interni + covarianza    
    if (opz.E_Pf_constraint) opt.m += 6*nLeg;

    if (opz.E_DV_cstr) opt.m += accumulate(nSeg_vector.begin(), nSeg_vector.end(), 0.) + nLeg; //ciascun DV è vincolato

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
    int nSeg;
    double tfin, Np1L;
    for (int iLeg = 0; iLeg < nLeg; iLeg++){
        nSeg = nSeg_vector[iLeg];
        tfin = tf_vector[iLeg];
        Np1L = accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg;
        for(int iSeg = 0; iSeg < nSeg + 1; iSeg++){
            int i0 = nVars*Np1L + iSeg*nVars;
            //int i0 = iLeg*(nSeg+1)*nVars + iSeg*nVars;

            if(iSeg < nSeg){ //DT
                opt.XL[i0] = 0.5*tfin/nSeg;
                opt.XU[i0] = tfin;
            }
            else{
                opt.XL[i0] = .0;
                opt.XU[i0] = .0;
            }

            for (int j = 0; j < 3; j++){ //pos
                opt.XL[i0 + 1 + j] = -2.0;
                opt.XU[i0 + 1 + j] = +2.0;
            }
            for (int j = 0; j < 3; j++){ //vel
                opt.XL[i0 + 1 + 3 + j] = -2.0;
                opt.XU[i0 + 1 + 3 + j] = +2.0;
            }
    
            for (int j = 0; j < 3; j++){ //DV
                if (opz.E_DV_cstr)
                {   
                    double LDV;
                    double UDV;
                    if ((opz.DV_RV_double == false) && iLeg > 0 && iSeg == 0){
                        LDV = UDV = .0;
                    }
                    else{
                        LDV = -opz.DVtot_single_max;
                        UDV = +opz.DVtot_single_max;
                    }
                    // opt.XL[i0+6+j] = -opz.DVtot_max/double(nSeg+1);
                    // opt.XU[i0+6+j] = +opz.DVtot_max/double(nSeg+1);
                    opt.XL[i0 + 1 + 6 + j] = LDV;
                    opt.XU[i0 + 1 + 6 + j] = UDV;
                }
                else
                {
                    opt.XL[i0 + 1 + 6 + j] = -1.0;
                    opt.XU[i0 + 1 + 6 + j] = +1.0;
                }
            }

            // KK
            if (opz.E_Pf_constraint){
                double LKK;
                double UKK;
                if ((opz.DV_RV_double == false) && iLeg > 0 && iSeg == 0){
                    LKK = UKK = .0;
                    }
                else{
                    LKK = -200.0;
                    UKK = +200.0;
                }
                //if(iLeg > 0) LKK = UKK = 0;
                for (int iRow = 0; iRow < 3; iRow++){       
                    for (int iCol = 0; iCol < 6; iCol++){
                        opt.XL[i0 + 1 + 6 + 3 + iRow*6 + iCol] = LKK;
                        opt.XU[i0 + 1 + 6 + 3 + iRow*6 + iCol] = UKK;
                    }
                }
            }
            else{
                for (int iRow = 0; iRow < 3; iRow++){       
                    for (int iCol = 0; iCol < 6; iCol++){
                        opt.XL[i0 + 1 + 6 + 3 + iRow*6 + iCol] = opt.XU[i0 + 6 + 3 + 1 + iRow*6 + iCol] = 0;
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
            r0_X = rRV;
            v0_X = vRV;
        }
        int i0 = nVars*Np1L;
        //int i0 = iLeg*(nSeg+1)*nVars;
        for (int j = 0; j < 3; j++){
            opt.XL[i0 + 1 + j] = opt.XU[i0 + 1 + j] = r0_X(j);
        }
        if((opz.v_RV_free == false)||(iLeg == 0)){                                 
            for (int j = 0; j < 3; j++){
            opt.XL[i0 + 1 + 3 + j] = opt.XU[i0 + 1 + 3 + j] = v0_X(j);
            }  
        }
        // else{
        //     for (int j=0; j<3; j++){
        //         opt.XL[i0 + 3 + j] = opt.XU[i0 + 3+j] = v0_X(j);
        //     } 
        // }

        // Vincoli: Residui + Condizioni al contorno
        // residuo r & v
        for (int iSeg = 0; iSeg < nSeg; iSeg++){
            int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + Np1L + iSeg*(nRes); // r, v, sigma_rf, sigma_vf, DV
            //int iG0 = iLeg*(nSeg*nRes + 6 + nSeg + 1) + iSeg*(nRes);
            if(opz.E_Pf_constraint) iG0 += 6*iLeg;

            //residuo_rv
            for (int j = 0; j < 6; j++){  
                opt.GL[iG0 + j] = opt.GU[iG0 + j] = 0; 
            }
        }

        // condizioni contorno finale
        int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + Np1L + nSeg*(nRes);
        //int iG0 = iLeg*(nSeg*nRes + 6 + nSeg + 1) + nSeg*(nRes);
        if (opz.E_Pf_constraint) iG0 += 6*iLeg;

        for (int i = 0; i < 6; i++){
            opt.GL[iG0 + i] = opt.GU[iG0 + i] = 0; 
        }

        if (opz.obj_func_switch == 1)   // vincola Pf e DV
        {
            // covarianza di posizione finale: Pr - Prdes < 0
            if (opz.E_Pf_constraint)
            {
                for (int j=0; j<6; j++)
                {
                    opt.GL[iG0 + 6 + j] = -par.Infty;   
                    opt.GU[iG0 + 6 + j] = 0;
                }
            } 
            // DVcstr: DV[k] - DV_max < 0
            if (opz.E_DV_cstr)
            {
                int nPf_cstr = 0;
                if (opz.E_Pf_constraint) nPf_cstr = 6;
                for (int j = 0; j < nSeg + 1; j++)
                {
                    opt.GL[iG0 + 6 + nPf_cstr + j] = -par.Infty; 
                    opt.GU[iG0 + 6 + nPf_cstr + j] = 0;
                }
            } 
        }
        else if (opz.obj_func_switch == 2)   // vincola DV e minimizza Pf
        {
            // DVcstr: DV[k] - DV_max < 0
            int nPf_cstr = 0;
            for (int j = 0; j < nSeg + 1; j++)
            {
                opt.GL[iG0 + 6 + nPf_cstr + j] = -par.Infty; 
                opt.GU[iG0 + 6 + nPf_cstr + j] = 0;
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
    
    int nSeg;
    double tfin, Np1L, dt;
    for (int iLeg = 0; iLeg < nLeg; iLeg++){
        nSeg = nSeg_vector[iLeg];
        vector<double> v_DVnorm(nSeg + 1); 
        vector<double> v_DVstd(nSeg + 1);

        //dt = (tf_vector[iLeg] - t0)/nSeg;

        Np1L = accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg;
        int iX = nVars*Np1L;

        dt = X[iX];
        r_km = Eigen::Map<Vector3d>(X + iX + 1, 3);  // Initial state is assigned by setting upper and lower boundaires in X equal to \tilde{x0} 
        v_km = Eigen::Map<Vector3d>(X + iX + 1 + 3, 3);       
        if (iLeg == 0){
            P_km = P0;  
        }        
        else{
            P_km = P_kp;
            if (opz.v_RV_free) v_km = v_kp;
        }                          
        DV_k = Eigen::Map<Vector3d>(X + iX + 1 + 6, 3);      
        KK_k = Eigen::Map<MatrixXd>(X + iX + 1 + 6 + 3, 3, 6); 
        KKaug_k << MatrixXd::Zero(3, 6), KK_k; 

        // apply DV
        v_kp = v_km + DV_k;
        P_kp = (eye6x6 + KKaug_k) * P_km * (eye6x6 + KKaug_k).transpose();
        covDV = KK_k*P_km*KK_k.transpose();

        if((opz.E_nav_std == 1)&&(iLeg == 0)){
            P_kp += KKaug_k*(RR_k*KKaug_k.transpose());
            covDV += KK_k*(RR_k*KK_k.transpose());
        }

        v_DVnorm[0] = DV_k.norm();  
        if (opz.DVstd_model == 0) 
            v_DVstd[0] = sqrt(covDV.trace());                         //sqrt(tr(DVstd))
        else    
            v_DVstd[0] = sqrt(covDV.eigenvalues().real().maxCoeff()); //max autovalore

        for (int iSeg = 0; iSeg < nSeg; iSeg++){
            // ottenuti per propagazione
            if (opz.cov_propagation_mode == 0)
                propagate_kepler_UT(r_km, v_kp, P_kp, dt, 1., Qd_k, r_k1_hat, v_k1m_hat, P_k1m_hat);
            else
                propagate_kepler_STM(r_km, v_kp, P_kp, dt, 1., Qd_k, r_k1_hat, v_k1m_hat, P_k1m_hat, opz.cov_propagation_nSub);

            int i0 = iX + (iSeg + 1)*nVars;
            dt = X[i0];
            Eigen::Map<Vector3d> r_k1m(X + i0 + 1, 3);  
            Eigen::Map<Vector3d> v_k1m(X + i0 + 1 + 3, 3);                                       
            P_k1m = P_k1m_hat;                                  // Hybrid: Covarianza propagata
            Eigen::Map<Vector3d> DV_k1(X + i0 + 1 + 6, 3);        
            Eigen::Map<MatrixXd> KK_k1(X + i0 + 1 + 6 + 3, 3, 6);        
            
            int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + Np1L + iSeg*(nRes);
            //int iG0 = iLeg*(nRes*nSeg + 6 + nSeg + 1) + nRes*(iSeg); 
            if (opz.E_Pf_constraint) iG0 += 6*iLeg;

            for (int j = 0; j < 3; j++)
                G[iG0 + j] = r_k1m[j] - r_k1_hat[j];
            for (int j = 0; j < 3; j++)
                G[iG0 + 3 + j] = v_k1m[j] - v_k1m_hat[j];

            // apply DV
            v_k1p = v_k1m + DV_k1;
            MatrixXd KKaug_k1(6, 6); KKaug_k1 << MatrixXd::Zero(3, 6), KK_k1;  
            MatrixXd P_k1p = (eye6x6 + KKaug_k1) * P_k1m_hat * (eye6x6 + KKaug_k1).transpose();  
            MatrixXd cov_DV_k1 = KK_k1*P_k1m*KK_k1.transpose();          

            // Navigation Errors
            if(opz.E_nav_std == 1)
                P_k1p += KKaug_k1*(RR_k*KKaug_k1.transpose());
                cov_DV_k1 += KK_k1*(RR_k*KK_k1.transpose());
            
            v_DVnorm[iSeg + 1] = DV_k1.norm(); 
            if (opz.DVstd_model == 0) 
                v_DVstd[iSeg + 1] = sqrt(cov_DV_k1.trace());                              // sqrt(tr(DVstd))
            else    
                v_DVstd[iSeg + 1] = sqrt(cov_DV_k1.eigenvalues().real().maxCoeff());      // max autovalore

            // update
            r_km = r_k1m; // note that: r_km = r_kp
            v_kp = v_k1p;
            P_kp = P_k1p;
        }

        // obj function
        double DVtot = accumulate(v_DVnorm.begin(), v_DVnorm.end(), 0.);
        double DVstd_tot = accumulate(v_DVstd.begin(), v_DVstd.end(), 0.);
        //     F = DVtot;
        F += DVtot + DVstd_tot * kstd;   

        // constraints
        //
        
        if (iLeg == nLeg - 1){
            r_G = opz.rf;
            v_G = opz.vf;
        }
        else{
            r_G = opz.rRV;
            v_G = opz.vRV;
        }
        
        // 1.1) Stato medio al tempo finale
        int iG0 = nRes*accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + 6*iLeg + Np1L + nSeg*(nRes);
        //int iG0 = iLeg*(nRes*nSeg + 6 + nSeg + 1) + nRes*nSeg;
        if (opz.E_Pf_constraint) iG0 += 6*iLeg;

        // if (iLeg == nLeg - 1){
        //     for (int j=0; j<3; j++)
        //         G[iG0 + j] = r_km[j] - r_G[j];

        //     for (int j=0; j<3; j++)
        //         G[iG0 + 3 +j] = v_kp[j] - v_G[j];
        // }
        // else{
        //     for (int j=0; j<3; j++)
        //         G[iG0 + j] = r_km[j] - r_G[j];

        //     for (int j=0; j<3; j++)
        //         G[iG0 + 3 +j] = 0;
        // }
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
            Pf_ev_r = P_kp.block(0,0,3,3).eigenvalues().real();
            Pf_ev_v = P_kp.block(3,3,3,3).eigenvalues().real();
            
            if (iLeg == nLeg - 1){
                sigma_r_G = opz.sigma_rf_des*prconv;
                sigma_v_G = opz.sigma_vf_des*prconv;
                }
            else{
                sigma_r_G = opz.sigma_rRV*prconv;
                sigma_v_G = opz.sigma_vRV*prconv;
            }
             
            for (int j = 0; j < 3; j++){
                //G[iG0 + 6 + j] = P_kp(j,j) - sigma_r_G; 
                G[iG0 + 6 + j] = opz.sf*(Pf_ev_r(j) - sigma_r_G);
                }
            for (int j = 3; j < 6; j++){
                //G[iG0 + 6 + j] = P_kp(j,j) - sigma_v_G; 
                G[iG0 + 6 + j] = opz.sf*(Pf_ev_v(j-3) - sigma_v_G); 
                }
        }

        // 1.3) DV_cstr
        if (opz.E_DV_cstr){
            int nPf_cstr = 0;
            if (opz.E_Pf_constraint) nPf_cstr = 6;

            for (int k = 0; k < nSeg + 1; k++){
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
    string fname_opz = opz_dir + "opz-rr-AFC.yaml"; 

    C_prb_RR_HMS prb(fname_opz);
    
    //--------------------------------//
    // inizializzazione
    //
    vector<int> nSeg_vector = prb.nSeg_vector;
    vector<double> tf_vector = prb.tf_vector;
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
    
    //--------------------------------//
    // Setup
    // 
    prb.setup_worph(Xguess);

    int MaxIter = 10000;
    // cout << "Jmax ?" << endl;
    // cin >> MaxIter;

    time_t tstart, tend; 

    if (MaxIter>0){
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
    ofstream trajMfile("../results/HMS-temp/trajMfile.dat");
    string fname_savings = "../results/HMS-temp/opt_X.dat";
    string out_opz("../results/HMS-temp/opz.yaml");
    double F, FLeg, dt, t, tL, tfin;
    int nSeg;
    F = 0;
    t = 0;

    MatrixXd P_kmf;
    P_kmf = MatrixXd::Zero(6,6);
    VectorXd r_k;
    VectorXd v_km;
    MatrixXd P_km;
    Vector3d r_k1, v_k1m;
    MatrixXd P_k1m(6, 6);
    MatrixXd Qd_k;
    MatrixXd RR_k;

    Vector3d v_kp;
    MatrixXd P_kp;
    MatrixXd cov_DV_i;

    vector<double> v_DVnorm, v_DVstd;

    C_KeplerMultiArc trajE, trajM;
    Vector3d r0_E, v0_E;
    r0_E << 43148032.1083303, 140976675.534384, -8649.81075546145;
    v0_E << -28.9752166461330, 8.61508360819854, 2.87602449144941e-06;
    r0_E /= prb.opz.rconv;
    v0_E /= prb.opz.vconv;
    trajE.add_keplerArc(C_KeplerArc(r0_E, v0_E, 0, 2.5*M_PI));
    trajE.print(trajEfile, 200);
    trajEfile.close();
    // trajM.add_keplerArc(C_KeplerArc(prb.rf_des, prb.vf_des, 0, 2.5*2.13*M_PI));
    // trajM.print(trajMfile, 200);
    // trajMfile.close();

    for(int iLeg = 0; iLeg < nLeg; iLeg++){
        tL = 0;
        nSeg = nSeg_vector[iLeg];
        tfin = tf_vector[iLeg];
        ofstream out_sommario(output_dir + "sommario" + to_string(iLeg + 1) + ".txt");
        string foldername = "../results/HMS-temp/fullsol" + to_string(iLeg + 1) + ".dat";
        ofstream out_DVs("../results/HMS-temp/DVs" + to_string(iLeg + 1) + ".dat");
        ofstream out_covEellipses("../results/HMS-temp/covEllipses" + to_string(iLeg + 1) + ".dat");
        ofstream out_sol(foldername);

        if (iLeg == 0){
            r_k = prb.r0;
            v_km = prb.v0;
            P_km = prb.P0;
        }
        else{
            r_k = prb.rRV;
            v_km = v_kp;
            P_km = P_kp;
        }
        Qd_k = prb.Qd_k;
        RR_k = prb.RR_k;

        C_KeplerMultiArc traj;
        traj.add_keplerArc(C_KeplerArc(r_k, v_km, 0, 0));

        for (int k = 0; k < nSeg + 1; k++){
            int i0 = nVars*(accumulate(nSeg_vector.begin(), nSeg_vector.begin() + iLeg, 0.) + iLeg) + k*prb.nVars;
            dt = prb.opt.X[i0];
            Eigen::Map<Vector3d> r_k(prb.opt.X + i0 + 1, 3);
            Eigen::Map<Vector3d> v_km(prb.opt.X + i0 + 1 + 3, 3);
            if ((prb.opz.v_RV_free)&&(k == 0)&&(iLeg > 0)) v_km = v_kp;
            Eigen::Map<Vector3d> DV_i(prb.opt.X + i0 + 1 + 6, 3);
            Eigen::Map<MatrixXd> KK_i(prb.opt.X + i0 + 1 + 6 + 3, 3, 6);
            MatrixXd KKaug_i(6, 6); KKaug_i << MatrixXd::Zero(3, 6), KK_i;

            // Apply DV
            v_kp = v_km + DV_i;
            P_kp = (eye6 + KKaug_i) * P_km * (eye6 + KKaug_i).transpose();
            cov_DV_i = KK_i*P_km*KK_i.transpose() * prb.prconv;   

            if((prb.opz.E_nav_std == 1)&&((k != 0)||(iLeg == 0))){
                P_kp += KKaug_i*(RR_k*KKaug_i.transpose());
                cov_DV_i += KK_i*(RR_k*KK_i.transpose());
            }

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

            out_DVs << fixed << setw(12) << setprecision(6) << t;
            for(int ii=0; ii<3; ii++) out_DVs << fixed <<setw(12) << setprecision(6) << r_k[ii];
            for(int ii=0; ii<3; ii++) out_DVs << fixed <<setw(12) << setprecision(6) << DV_i[ii]; 
            out_DVs << endl;

            print_covariance(r_k, v_kp, P_km * prb.prconv, out_covEellipses);
            print_covariance(r_k, v_kp, P_kp * prb.prconv, out_covEellipses);

            out_sommario << fixed << setw(2) << setprecision(8) << k << scientific << setw(2) << setprecision(5) << P_km.diagonal().transpose() << endl;
            out_sommario << fixed << setw(2) << setprecision(8) << k << scientific << setw(2) << setprecision(5) << P_kp.diagonal().transpose() << endl;
            out_sommario << endl;

            if (k != nSeg){
                // propagate
                
                // Vector3d r_k1_hat, v_k1m_hat;  MatrixXd P_k1m_hat;
                if (prb.opz.cov_propagation_mode == 0)
                    propagate_kepler_UT(r_k, v_kp, P_kp, dt, 1., Qd_k, r_k1, v_k1m, P_k1m);
                else
                    propagate_kepler_STM(r_k, v_kp, P_kp, dt, 1., Qd_k, r_k1, v_k1m, P_k1m, prb.opz.cov_propagation_nSub);

                traj.add_keplerArc(C_KeplerArc(r_k, v_kp, t, dt));
            }
            else{
                r_k1 = r_k;
                v_k1m = v_kp;
                P_k1m = P_kp;
                P_kmf = P_km;
                traj.add_keplerArc(C_KeplerArc(r_k, v_kp, t, 0));
            }

            // update
            P_km = P_k1m;
            t += dt;
            tL += dt;
        }

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
        out_sommario << "DTLeg = " << tL << endl;
        out_sommario << "DT = " << tfin << endl;

        v_DVnorm.clear(); v_DVstd.clear();
        traj.print(out_sol, 100);
        if (iLeg == nLeg - 1){
            out_sommario << endl;  
            out_sommario << "F = " << F << endl;
        }

        out_sommario.close();
        out_sol.close();
        out_DVs.close();
        out_covEellipses.close();
    }

    // system("cd ..; gnuplot plot-files/plot_traj_RR_HMS.plt");

    // Save YP to Disk
    prb.save_sol(fname_savings);
    prb.opz.emit(out_opz);
    
    system("cd ..; gnuplot > /dev/null 2>&1 plot-files/plot_traj_RR_HMS.plt");

    string ans = "S";
    //cout << "Copy solution in dbSOL? (S/n)" << endl;
    //cin >> ans;
    if (ans.compare("S") == 0 || ans.compare("s") == 0){
        std::ifstream if_a("../results/HMS-temp/fullsol1.dat", std::ios_base::binary);
        std::ifstream if_b("../results/HMS-temp/fullsol2.dat", std::ios_base::binary);
        std::ofstream of_c("../results/HMS-temp/fullsol.dat", std::ios_base::binary);

        of_c << if_a.rdbuf() << if_b.rdbuf();

        string cmd;
        cmd = "mkdir -p ../dbSOL/" + prb.opz.output_folder; system(cmd.c_str()); //cout << cmd << endl;
        cmd = string("cp -RT ") + "../results/HMS-temp" + " ../dbSOL/" + prb.opz.output_folder;  system(cmd.c_str()); //cout << cmd << endl;
        // cmd = "cp -RT " + output_dir + " ../dbSOL/"+prb.opz.output_folder;  system(cmd.c_str()); cout << cmd << endl;
    }
    cout << "It took " << difftime(tend, tstart) << " second(s) to solve the ROCP."<< endl;
}