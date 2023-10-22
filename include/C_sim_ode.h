#ifndef H_C_SIM_ODE
#define H_C_SIM_ODE

# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <fstream>
# include <iostream>
#include <vector>

#include<eigen3/Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

/************************************************************************************
*
*
*										Core
*			FUNC f, //void f ( double t, double y[], double yp[] ),
*/
template<typename FUNC>
int ode_driver_eigen(FUNC odesys, const VectorXd& y0, VectorXd& T, MatrixXd& Y, const double& abserr, const double& relerr)
{
 
 
	int iflag;
	int iwork[5];
	//int neqn = 2;
	const int neqn = y0.size();
	//int step_num = 12;
	int NCOL = T.size();


	double t;
	double tout;
	double *work;
	


 

	iflag = -1;   //normale � +1


	
	// Set first point 
	//t = 0.;
	//y[0] = 1.0;
	//y[1] = 0.0;
	t = T(0);
	//y = new double[neqn];
	//y = y0.data();
	
	//double *y;
	//Eigen::Map<VectorXd>(y, y0.rows()) = y0;
	double* y = new double[neqn];
	for (int i = 0; i < neqn; i++)	y[i] = y0(i);


	//** Save first point (y0) 	
	//cout << "  " << setw(8) << t 	<< "  " << setw(14) << y[0]		<< "  " << setw(14) << y[1] << "\n";
	Y.col(0) = y0;

	work = new double[100 + 21 * neqn];

	for (int j = 1; j < NCOL; j++)
	{
		tout = T(j);

		ode(odesys, neqn, y, t, tout, relerr, abserr, iflag, work, iwork); //iflag == 2 -> correct 


		if (iflag != 2)
		{
			cout << "\n";
			cout << "TEST01 - Fatal error!\n";
			cout << "  ODE returned IFLAG = " << iflag << "\n";
			return 1;
		}
		
		// non superare l'ultimo passo
		//if (j == NCOL-2) 
		iflag = -2;
		
		// save
		//cout << "  " << setw(8) << t 	<< "  " << setw(14) << y[0]		<< "  " << setw(14) << y[1] << "\n";
		Y.col(j) = Eigen::Map<VectorXd>(y, neqn);


	}
	delete[] work;
	delete[] y;
	return 0;
}



class C_simulazione_ODE
{
public:
	int NEQ;
	int NCOL;
	int ICMP;
	VectorXd y0;		//c.i.
	VectorXd T;			//output - Tempo
	MatrixXd Y;			//output - Traiettoria  NEQ x NCOL
	double ATOL, RTOL;	//tolleranze di integrazione


	C_simulazione_ODE(int NEQ_, int NCOL_, int ICMP_, const double& ATOL_, const double& RTOL_) : NEQ(NEQ_), NCOL(NCOL_), ATOL(ATOL_), RTOL(RTOL_), ICMP(ICMP_)
	{
		y0 = VectorXd(NEQ);
		T = VectorXd(NCOL);
		Y = MatrixXd(NEQ, NCOL);


		//Tspan - lineare
		if (ICMP > 0)
			T.setLinSpaced(NCOL, ICMP - 1, ICMP);
		else
		{
			std::cout << "errore setLinSpaced" << std::endl;
			T.setLinSpaced(NCOL, 0., 1.);
		}

	}




	void SetY0_Sim(VectorXd y0_in)
	{
		if (y0_in.size() != NEQ)
		{
			std::cout << "errore nelle dimensioni" << std::endl
				<< "y0.size()  = " << y0.size() << std::endl
				<< "y0_in.size() = " << y0_in.size() << std::endl;
		}
		else
			y0 = y0_in;
	};


	template<typename FUNC>
	int Start_Sim(FUNC odesys)
	{
		int ierr = ode_driver_eigen(odesys, y0, T, Y, ATOL, RTOL);		//non resituisce mai problemi 		


		if (ierr > 0)
		{
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
		return ierr; //usicita corretta
	};




	void Print_Sim(std::ostream& os)
	{
		/*TYPE(Simulazione) ::sim
		LOGICAL, intent(in), OPTIONAL::first, last

		IF((.NOT.present(first)).OR.(first)) THEN
		OPEN(UNIT = 10, FILE = 'example1.esp')
		WRITE(10, *) "/td 'xyy'"
		ENDIF*/
		for (int JJ = 0; JJ < NCOL; JJ++)
		{
			os << std::fixed << std::setw(10) << std::setprecision(6) << T(JJ) << "\t" << Y.col(JJ).transpose() << std::endl;
		}
	};


};


































#endif //H_C_SIM_ODE