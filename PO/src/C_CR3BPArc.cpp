#include "../include/C_CR3BPArc.h"

void C_CR3BPMultiArc::add_CR3BPArc(C_CR3BPArc arc){
    v_CR3BPArcs.push_back(arc);
}

void C_CR3BPMultiArc::print(ostream& os, int npt){
    int nArcs = v_CR3BPArcs.size();
    for (int iArc = 0; iArc < nArcs; iArc ++){
        if (iArc == 0)
            v_CR3BPArcs[iArc].print(os, npt, true);
        else
            v_CR3BPArcs[iArc].print(os, npt);
        os << endl;
    }
}

void C_CR3BPArc::get_final_state(Vector3d &rf, Vector3d &vf)
{
    if (not(final_state)){
        if (fabs(tf - t0) < 1e-8){
            _rf = r0;
            _vf = v0;
        }
        else{
            VectorXd s0(6), Yf(6); s0 << r0, v0;
            C_simulazione_ODE sim(6, 2000, 1, 1e-12, 1e-14);
            sim.SetY0_Sim(s0);
            EoM_CR3BP_ODE EoM(mu, (tf - t0));
            sim.Start_Sim(EoM);
            Yf = sim.Y.rightCols(1);
            _rf = Yf.segment(0,3);
            _vf = Yf.segment(3,3);
        }
        final_state = true;
    }
    rf = _rf;
    vf = _vf;
};

void C_CR3BPArc::print(ostream &os, int npt, bool header_on){

    //const int npt = 10;
    int w = 20;
    double tconv = 1; //5.03626e+06;
    double rconv = 1; //1.49868e+08;
    double vconv = 1; //29.7578;

    if (false)
        os << endl;
        os << fixed << std::setw(w) << setprecision(6) << t0*tconv
           << fixed << std::setw(w) << setprecision(6) << r0(0)*rconv
           << fixed << std::setw(w) << setprecision(6) << r0(1)*rconv
           << fixed << std::setw(w) << setprecision(6) << r0(2)*rconv
           << fixed << std::setw(w) << setprecision(6) << v0(0)*vconv
           << fixed << std::setw(w) << setprecision(6) << v0(1)*vconv
           << fixed << std::setw(w) << setprecision(6) << v0(2)*vconv
           << endl;

    if (abs(tf - t0) > 1e-8){
        double dt_inner = (tf - t0) / double(npt - 1);
        C_simulazione_ODE sim(6, 2000, 1, 1e-12, 1e-14);
        Vector3d r_k, v_k;
        VectorXd s0(6), Yf(6); s0 << r0, v0;
        sim.SetY0_Sim(s0);
        for (int ipt = 1; ipt < npt; ipt++){
            double t_k = t0 + ipt * dt_inner;

            EoM_CR3BP_ODE EoM(mu, t_k - t0);
            sim.Start_Sim(EoM);
            Yf = sim.Y.rightCols(1);
            r_k = Yf.segment(0,3);
            v_k = Yf.segment(3,3);

            // save
            os << fixed << std::setw(w) << setprecision(6) << t_k*tconv
               << fixed << std::setw(w) << setprecision(6) << r_k(0)*rconv
               << fixed << std::setw(w) << setprecision(6) << r_k(1)*rconv
               << fixed << std::setw(w) << setprecision(6) << r_k(2)*rconv
               << fixed << std::setw(w) << setprecision(6) << v_k(0)*vconv
               << fixed << std::setw(w) << setprecision(6) << v_k(1)*vconv
               << fixed << std::setw(w) << setprecision(6) << v_k(2)*vconv
               << endl;
        }
    }
};