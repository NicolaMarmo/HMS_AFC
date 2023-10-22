#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include "myRandom.h"
#include <cmath>

using namespace std;


// v1.1
//  ++ add :  genotypeDistance()
// 
//  ++ add : fullDisplayDISK  (save a tmp file with Gbest @ full precision)
//
// v1.2	@ 16/10/2017
//  ++ add :  comparison by alpha level (satisfaction level) 



 

class C_agent
{
public:
	//--------------//
	vector<double> x;		 /*!< genotipe */
	double fit;				 /*!< fitness */
	vector<double> cstrViol; /*!< constraint violation */
	double weighCstrViol;	 /*!< total weighted constraint violation */
	double fitAug;			 /*!< fit + penality */ 
	//--------------//


	C_agent() {} ;
	C_agent(const int& d) { x.resize(d); };
	C_agent(vector<double> x_in) : x(x_in) {};

	//member functions
	//void C_agent::fullDisplayDISK(ostream &os) const;
	void fullDisplayDISK(ostream &os) const;


	//friend bool operator<(C_agent a_left, C_agent a_right);
	friend C_agent rndGenAgent(vector<double> lb, vector<double> ub);

	//cmp operator
	friend bool cmp_byPenality(C_agent a_left, C_agent a_right);
	friend bool cmp_byRule1(C_agent a_left, C_agent a_right);
	friend bool cmp_by_weightCstrViol(C_agent a_left, C_agent a_right);
	friend double genotypeDistance(const C_agent& a_left, const C_agent& a_right);

private:

};





// --------------------------------- //
//			Classe Jagent			 //
// --------------------------------- //

class C_jagent : public C_agent{
public:
	double F; // = 0.5;			/*!<	Mutation parameter */
	double Cr; // = 0.8;		/*!<	Mutation parameter */
	int iStrategy; // = 1;		/*!<	Mutation parameter */


	//Costructors
	C_jagent();
	C_jagent(C_agent a_in);		//inverse-cast
	
	//friend functions
	void fullDisplay(ostream &out) const;
	//void fullDisplayDISK(ostream &out) const;
	friend void assignAgent(C_agent a_in, C_jagent& aj_out);


private:
};


bool cmp_byPenality(C_agent a_left, C_agent a_right);
bool cmp_by_weightCstrViol(C_agent a_left, C_agent a_right);
bool cmp_byRule1(C_agent a_left, C_agent a_right);


struct t_cmp_epsilon
{
public:
	double epsLvL;


	t_cmp_epsilon(const double& epsLvL) : epsLvL(epsLvL) {};

	bool operator()(C_agent a_left, C_agent a_right)
	{
		/*!
		*	comparison function object (i.e. an object that satisfies the requirements of Compare)
		*	which returns ​true if the first argument is less than the second.
		*
		*	Regola:
		*	1) Fra due individui feasible -> quello con la fitness migliore
		*	2) Fra un feasible and a non-feasible -> quello feasible
		*	3) Fra due non-feasible, quello con la minore totalWeightedCstrViolation
		*
		*		totalWeightedCstrViolation = cstrViol[i]/(eps+cstrViolMax_G[i])
		*
		*		cstrViol[i] = violazione del vincolo i-th
		*		cstrViolMax_G[i] = violazione massima del vincolo i-th alla generazione precedente (G)
		*
		*	Return:
		*		true = xi < ui,		false = xi > ui
		*/


		if (a_left.weighCstrViol <= epsLvL && a_right.weighCstrViol <= epsLvL)
		{
			//caso 1
			if (a_left.fit < a_right.fit)
				return true; //select xi
			else
				return false; //select ui
		}
		else if (a_left.weighCstrViol > epsLvL && a_right.weighCstrViol > epsLvL)
		{
			//caso 3
			if (a_left.weighCstrViol < a_right.weighCstrViol)
				return true; //select xi
			else
				return false; //select ui
		}
		else //solo uno dei due è feasible, l'altro no
		{
			//caso 2:
			if (a_left.weighCstrViol < epsLvL)
				return true; //select xi
			else
				return false; //select ui
		}


	};

};
 