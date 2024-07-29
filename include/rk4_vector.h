#ifndef H_RK4_EIGEN
#define H_RK4_EIGEN

/*
ToDO:	il metodo di integrazione deve resituire un codice se � tutto ok o si � silurato qualcosa.
*/


# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <fstream>
# include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace std;

/**
  * 
  *   richide di scrivere il sistema di equazioni ODE come una funzione/FUNTORE con signature
  *       VectorXd   Func(double, VectorXd)
  *
  *	
  */

//VectorXd myfunz(double t, VectorXd y);


// assume esista un funtore "Func = double *f(double t, int m, double u[])"
template<typename Func>
VectorXd rk4_eigvec(Func f, double t0, VectorXd u0, double dt)
{

	///*   VectorXd   Func(double, VectorXd)

	//int m = u0.size();
	//VectorXd f0(m), f1(m), f2(m), f3(m);	
	//VectorXd u(m), u1(m), u2(m), u3(m);


	//
	//  Get four sample values of the derivative.
	//
	VectorXd f0 = f(t0, u0);

	double t1 = t0 + dt / 2.0;
	VectorXd u1 = u0 + dt * f0 / 2.0;
	VectorXd f1 = f(t1, u1);

	double t2 = t0 + dt / 2.0;
	VectorXd u2 = u0 + dt * f1 / 2.0;
	VectorXd f2 = f(t2, u2);

	double t3 = t0 + dt;	
	VectorXd u3 = u0 + dt * f2;
	VectorXd f3 = f(t3, u3);
	//
	//  Combine them to estimate the solution.
	//
	VectorXd u = u0 + dt * (f0 + 2.0 * f1 + 2.0 * f2 + f3) / 6.0;

	//}
	////
	////  Free memory.
	////
	//delete[] f0;
	//delete[] f1;
	//delete[] f2;
	//delete[] f3;
	//delete[] u1;
	//delete[] u2;
	//delete[] u3;

	return u;
};



template<typename FUNC>
void rk4_driver_eigen(FUNC odesys, const VectorXd& y0, double t0, double  tf, int nStep,	//in
	VectorXd& T, MatrixXd& Y)			//out
{
	// rk4_driver ( odesys, y0, t0, tf, nStep, T, Y)
	/*

	y0 [NY]
	t0
	tf	double
	T [nStep+1]
	Y [NY][nStep+1]
	*/

	int M = nStep + 1;
	int NY = y0.size();

	T.setLinSpaced(nStep+1, t0, tf);
	Y.setZero(NY, M);


	double dt = (tf - t0) / double(nStep);

	// Save @ starting time 
	T(0) = t0;
	Y.col(0) = y0;




	for (int istep = 0; istep < nStep; istep++)
	// integro un passo con RK4
	{	
		// I.C.	
		VectorXd y = Y.col(istep);
		double t = T(istep);		
		double dt = T(istep + 1) - T(istep);

		// stepper
		VectorXd y1 = rk4_eigvec(odesys, t, y, dt);

		//update
		Y.col(istep + 1) = y1;
	}
};



/******************************************************************************/

template<typename FUNC>
int rk4_driver_eigen(FUNC odesys, const VectorXd& y0, const VectorXd& T, //in
	MatrixXd& Y)			//out
{


	int ierr = 0; //         0 = OK, 1= ERRORE

	// rk4_driver_V2 ( odesys, y0, t0, tf, nStep, T, Y)
	/*

	y0 [NY]
	t0
	tf	double
	T [nStep+1]
	Y [NY][nStep+1]
	*/
	int M = T.size();
	int nStep = M - 1;	
	int NY = y0.size();

	
	
	
	Y.setZero(NY, M);	// Inizialize
	Y.col(0) = y0;		// Save @ starting time 

	
//	double dt = (tf - t0) / double(nStep);
	for (int istep = 0; istep < nStep; istep++)
		// integro un passo con RK4
	{
		// I.C.	
		VectorXd y = Y.col(istep);
		double t = T(istep);
		double dt = T(istep + 1) - T(istep);

		// stepper
		VectorXd y1 = rk4_eigvec(odesys, t, y, dt);

		//update
		Y.col(istep + 1) = y1;
	}

	for (int i = 0; i < Y.rows(); i++)
	{
		for (int j = 0; j < Y.cols(); j++)
		{
			if (isnan(Y(i, j)))
			{
				ierr = 1;
				// cout << "c'� un NAN" << endl;
			}
		}
	}
	/*if (Y.any().isnan())
		cout << "c'� un NAN" << endl;*/


	return  ierr;

};



/******************************************************************************/


class C_simulazione
{
public:
	int NEQ;
	int NCOL;
	int ICMP;					
	VectorXd y0;		//c.i.
	VectorXd T;			//output - Tempo
	MatrixXd Y;			//output - Traiettoria  NEQ x NCOL
	double ATOL, RTOL;	//tolleranze di integrazione

	

	C_simulazione(int NEQ_, int NCOL_, int ICMP_, const double& ATOL_, const double& RTOL_) : NEQ(NEQ_), NCOL(NCOL_), ATOL(ATOL_), RTOL(RTOL_), ICMP(ICMP_)
	{
		y0 = VectorXd(NEQ);
		T = VectorXd(NCOL);
		Y = MatrixXd(NEQ, NCOL);
		
		//Tspan - lineare
		if (ICMP > 0)
			T.setLinSpaced(NCOL, ICMP - 1, ICMP);
		else
		{
			cout << "errore setLinSpaced" << endl;
			T.setLinSpaced(NCOL, 0., 1.);
		}
			
	}




	void SetY0_Sim(VectorXd y0_in)
	{
		if (y0_in.size() != NEQ)
		{
			std::cout << "errore nelle dimensioni" << endl
				<< "y0.size()  = " << y0.size() << endl
				<< "y0_in.size() = " << y0_in.size() << endl;
		}
		else		
			y0 = y0_in;
	};


	template<typename FUNC>
	int  Start_Sim(FUNC odesys)
	{
		//IMPLICIT REAL * 8 (a - h, o - z)
		//!arg
		//EXTERNAL F
		//TYPE(Simulazione) ::sim
		//REAL * 8 ::YP(:)
		//!temp Array
		//REAL * 8, allocatable::Y_(:)
		//!temp Scalar
		//REAL * 8 ::T_, TOUT
		//!---------- -

		//	!da fare :
		//!1) se sim non � stata inizializzata->segnala errore
		//	!2) se iflag.ne. 2->segnala errore


		//call integratore(sim%NEQ, F, sim%T, sim%Y0, YP, sim%ATol, sim%RTol, sim%Y, sim%icmp, iflag)
		int ierr = rk4_driver_eigen(odesys, y0, T, Y);		//non resituisce mai problemi 

		if (ierr > 0)
		{
			/*
			cout << "errore nell'integrazione" << endl
				<< "subroutine ' Start_Sim '" << endl
				<< "metto un numero alto a caso per far riprovare il conto" << endl;
				*/
			Y.fill(1000.);
			return 1;
		}
		/*if (iflag.ne.2) then
			WRITE(6, *) "errore nell'integrazione"
			WRITE(6, *) "subroutine ' Start_Sim '"
			WRITE(6, *) "metto un numero alto a caso per far riprovare il conto"
			sim%Y = 1000d0
			!pause
			ENDIF
			*/
		return ierr; //uscita corretta
	};

	void Print_Sim(ostream& os)
	{
			/*TYPE(Simulazione) ::sim
			LOGICAL, intent(in), OPTIONAL::first, last

			IF((.NOT.present(first)).OR.(first)) THEN
			OPEN(UNIT = 10, FILE = 'example1.esp')
			WRITE(10, *) "/td 'xyy'"
			ENDIF*/
		for (int JJ = 0; JJ < NCOL; JJ++)
		{
			os << fixed << setw(10) << setprecision(6) << T(JJ) << "\t" << Y.col(JJ).transpose() << endl;
		}					
	};


};





/****************************************************************/
/*							Test								*/
/*																*/
/****************************************************************/












#endif 