#include "../include/aux_covariance.h"

void print_covariance(Vector3d r, Vector3d v, MatrixXd P, ostream &os)
{
    const double chi_sq = 1; //5.991; // 95% ellipse
        Eigen::EigenSolver<MatrixXd> es(P.block(0,0,2,2));

            MatrixXd A = es.eigenvectors().real();
            A.col(0) = A.col(0)*sqrt(chi_sq *es.eigenvalues()(0).real());
            A.col(1) = A.col(1)*sqrt(chi_sq *es.eigenvalues()(1).real());

            int npt=21;
            VectorXd vth = VectorXd::LinSpaced(npt,0, 2.*M_PI);
            MatrixXd B(2,npt); 
            B.row(0) = vth.array().sin();
            B.row(1) = vth.array().cos();
            MatrixXd CC = A*B;

            for (int ipt=0; ipt<npt; ipt++)
            {
                os //<< "# iSeg = " << k << endl
                    << fixed <<setw(12) << setprecision(6) << r[0] <<"\t"
                    << fixed <<setw(12) << setprecision(6) << r[1] <<"\t"
                    << scientific <<setw(12) << setprecision(6) << CC(0, ipt) <<"\t"
                    << scientific <<setw(12) << setprecision(6) << CC(1, ipt) <<"\t"
                    << endl;
            }
            os << endl << endl;

}

void print_DVs(double t, Vector3d r, Vector3d DV, ostream &os)
{
    os << fixed <<setw(12) << setprecision(6) << t;
    for(int ii=0; ii<3; ii++) os << fixed <<setw(12) << setprecision(6) << r[ii];
    for(int ii=0; ii<3; ii++) os << fixed <<setw(12) << setprecision(6) << DV[ii]; 
    os << endl;

}
     

