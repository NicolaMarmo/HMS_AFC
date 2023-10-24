#include "../include/C_KeplerArc.h"

C_KeplerMultiArc::C_KeplerMultiArc(Vector3d _r0, Vector3d _v0, VectorXd v_t, MatrixXd DVs){
    // asssume DV \in [3 x N+1]
    // T in [N+1]
    
    Vector3d r0 = _r0;
    Vector3d v0 = _v0;

    int nSeg = v_t.size() - 1;

    // primo arco ha durata nulla
    add_keplerArc(C_KeplerArc(r0,v0, v_t[0], 0.));

    double t0;
    for (int iSeg = 0; iSeg< nSeg; iSeg++){
        double tEnd = v_t[iSeg+1];
        Vector3d v0p = v0 + DVs.col(iSeg);
        this->add_keplerArc(C_KeplerArc(r0,v0p,t0, tEnd-t0));
        (this->v_keplerArcs.end()-1)->get_final_state(r0,v0);
        t0 = tEnd;
    }


    // utlimo arco ha durata nulla
    Vector3d v0p = v0 + DVs.col(nSeg);
    add_keplerArc(C_KeplerArc(r0,v0p, v_t[nSeg], 0.));
}

void C_KeplerMultiArc::add_keplerArc(C_KeplerArc arc){
    v_keplerArcs.push_back(arc);
}

void C_KeplerMultiArc::print(ostream& os, int npt){
    int nArcs = v_keplerArcs.size();
    for (int iArc=0; iArc<nArcs; iArc ++)
    {
        if (iArc == 0)
            v_keplerArcs[iArc].print(os, npt, true);
        else
            v_keplerArcs[iArc].print(os, npt);
        //os << endl;
    }
    
}

void C_KeplerArc::get_final_state(Vector3d &rf, Vector3d &vf){
    if (not(final_state))
    {
        if (fabs(tf - t0) < 1e-8)
        {
            _rf = r0;
            _vf = v0;
        }
        else
        {
            propagateKEP_U(r0, v0, tf - t0, 1., _rf, _vf);
        }
        final_state = true;
    }
    rf = _rf;
    vf = _vf;
};

void C_KeplerArc::print(ostream &os, int npt, bool header_on){
    //const int npt = 10;

    if (header_on)
        os << endl;
        os << fixed << std::setw(12) << setprecision(6) << t0
           << fixed << std::setw(12) << setprecision(6) << r0(0)
           << fixed << std::setw(12) << setprecision(6) << r0(1)
           << fixed << std::setw(12) << setprecision(6) << r0(2)
           << fixed << std::setw(12) << setprecision(6) << v0(0)
           << fixed << std::setw(12) << setprecision(6) << v0(1)
           << fixed << std::setw(12) << setprecision(6) << v0(2)
           << endl;

    if (abs(tf - t0) > 1e-8)
    {
        double dt_inner = (tf - t0) / double(npt - 1);
        for (int ipt = 1; ipt < npt; ipt++)
        {
            double t_k = t0 + ipt * dt_inner;

            Vector3d r_k, v_k;
            propagateKEP_U(r0, v0, t_k - t0, 1., r_k, v_k);

            // save
            os << fixed << std::setw(12) << setprecision(6) << t_k
               << fixed << std::setw(12) << setprecision(6) << r_k(0)
               << fixed << std::setw(12) << setprecision(6) << r_k(1)
               << fixed << std::setw(12) << setprecision(6) << r_k(2)
               << fixed << std::setw(12) << setprecision(6) << v_k(0)
               << fixed << std::setw(12) << setprecision(6) << v_k(1)
               << fixed << std::setw(12) << setprecision(6) << v_k(2)
               << endl;
        }
    }
};

