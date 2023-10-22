#pragma once
/*!

versione 1.3
requires: jeval.h v1.3+

v1.2
	++ add  - enforceBoxCstr: now support circular box constraints
	++ add  - pruneByClustering
	++ add  - OPENMP support:
	use "#define OMP_ON"  to enable/disable omp
	++ add  - mutation strategy GA_like_1 e GA_like_2
	++ modified - DE algorithm default parameters are set here

v1.3
	++ add - epidemia (or internal reset with elitism)

v 1.4
	++ add - popReduction

*/

#include <omp.h>

#include <assert.h>
#include <vector>
#include <cmath>

#include <algorithm>
#include "permutazioni.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "c_agent.h"
#include "jeval.h"

using namespace std;
 



//omp parameters
#define OMP_ON 1  //1 = active; 0 = inactive
#define N_ACTIVE_THREADS 32 //8



// DE algorithm parameters
const int nStrategyMax = 7;
const int selectionType_DEFAULT = 0;	// 0 = penality method
										// 1 = Rule1	(feasible-first method)

enum Enum_selectTionType
{
	penalty,
	feasibleRule,
	EPSILON_RULE
};

// **************************************************************************************** //
//				parametri modificabili x algoritmo Differential Evolution					//


//-	constraint handling

const int weighCstrMode = 0;			// 0 = non weighted; 1 = weighted


// -	 adaptation modes
const int adaptMode[] = { 2, 2, 1 }; //{ 2, 2, -4 };	// [F e Cr, strategy] adaptation modes
										// **  F **
										// 0 = noAdapt | 1 = popAdap | 2 - selfAdapt
										//
										// **  CR  **
										// 0 = noAdapt | 1 = popAdap | 2 - selfAdapt
										// 
										// **  S  **
										// [1:nStrategyMax]->quella strateguia costante
										//  0 = init a caso, costante(const across pop)
										// -1 = init a caso, circolare(const across pop)
										// -2 = init a caso, popAdap : rnd(const across pop)
										// -3 = init a caso, popAdap : unifRoulette(const across pop)
										// -4 = init a caso, selfAdapt : rnd(variabile across pop)
										// -5 = init a caso, selfAdapt : unifRoulette(variabile across pop)

const bool isPruningOn = 0;				// 0 = disable pruning | 1 = enable
const int nGenPruning = 500; //ogni quante generazioni fare il pruning

// -	terminal conditions
const int nGen_MIN_DEFAULT = 10000;		//	default max-generation limit
const int nGen_MAX_DEFAULT = 50000;		//	default max-generation limit
const double x_Gbest_tol = 1e-7;
const double fit_Gbest_tol = 1e-5;

const int inactiveGen_x_GbestMAX = 12000;
const int inactiveGen_fit_GbestMAX = 12000;
const int Max_FES = 500000 * 100;


const int nElite = 3;
const int nEpidemiaMax = 0;
const int nGenEpidemia = 500;
const double percEpidemia = 1.;		//perc di pop colpita da epidemia [0, 1.]

// **************************************************************************************** //

//Display parameters
extern ofstream hof_log;		//hof = hall of fame
extern ofstream tmpSol;
const int nGen_hof_Flush = 50;	//ogni quanto effettuare il flush dei dati hof


extern int FES;


ostream& operator<<(ostream& os, const C_agent& a);

//bool operator<(C_agent a_left, C_agent a_right);
bool cmp_byPenality(C_agent a_left, C_agent a_right);
bool cmp_byRule1(C_agent a_left, C_agent a_right);



void eval_Fit_and_Cstr(C_agent& agente);
void eval_weighCstrViol(C_agent& a_ui, const vector<double>& maxCstrViol);
//void eval_fitAug(C_agent& a_ui, const vector<double>& maxCstrViol, const int& nGen, const double& cstrTol = 1e-6);
void eval_fitAug(C_agent& a_ui, const vector<double>& maxCstrViol, const int& nGen);




C_agent engine_de(int d, int nCstr, vector<double> lb, vector<double> ub);

C_agent selectAgent(const C_agent& a_xi, const C_agent& a_ui, vector<double> cstrViolMax_G = vector<double>(1, -1.));






C_agent engine_de_OO(MinCstrPblmbAbsP ptr_Problem);

C_jagent engine_jde_OO(MinCstrPblmbAbsP ptr_Problem, int _np = -1, int _nGenMAX = -1, int _selectionType = -1);



template <typename T_agent>
T_agent findBestOfDeme(vector<T_agent> pop, Enum_selectTionType eSelectionType, const double epsLvL = 1e-6)
{
	T_agent Gbest;
	//	cerca il migliore del gruppo, according to a comparison operator

	switch (eSelectionType)
	{
	case feasibleRule:
		//cout << "cmp_byRule1" << endl;
		Gbest = *min_element(pop.begin(), pop.end(), cmp_byRule1);
		break;
	case penalty:
		//cout << "penalty" << endl;
		Gbest = *min_element(pop.begin(), pop.end(), cmp_byPenality);
		break;
	case EPSILON_RULE:
		//cout << "t_cmp_epsilon" << endl;
		Gbest = *min_element(pop.begin(), pop.end(), t_cmp_epsilon(epsLvL));
		break;
	}
	return Gbest;
};

template <typename T_agent>
T_agent findBestOfDeme(vector<T_agent> pop, int _selectionType, const double epsLvL = 1e-6)
{
	// epsLvL opzionale, solo per comparazione tipo epsilon


	T_agent Gbest;
	//	cerca il migliore del gruppo, according to a comparison operator
	if (_selectionType == 0)
		Gbest = *min_element(pop.begin(), pop.end(), cmp_byPenality); 
	else if(_selectionType == 1)
		Gbest = *min_element(pop.begin(), pop.end(), cmp_byRule1);	
	else if (_selectionType == 2)
	{
		//double epsLvL = 1e-6;
		Gbest = *min_element(pop.begin(), pop.end(), t_cmp_epsilon(epsLvL));
	}
		

	return Gbest;
	//C_agent findBestOfDeme(vector<C_agent> pop, int _selectionType)
};


template <typename T_agent>
double eval_epsLvL0(vector<T_agent> pop)
{
	sort(pop.begin(), pop.end(), cmp_by_weightCstrViol);   // sortDemeBy_weightCstrViol
	
	int N0 = ceil(0.05*pop.size());
	return pop[N0].weighCstrViol;
};




template<typename T_agent>
vector<double> getMaxCstrViol(const vector<T_agent>&pop){
	/// T_agent � un template che assume T = C_agent
	///  ma va bene anche per C_jagent
	///
	// 

	int nCstr = pop[0].cstrViol.size();
	int np = pop.size();

	vector<double> maxCstrViol(nCstr, 0);

	for (int j = 0; j < pop[0].cstrViol.size(); j++)
	{
		vector<double> vCstr_j(np, 0);
		for (int i = 0; i < np; i++)
			vCstr_j[i] = pop[i].cstrViol[j];

		maxCstrViol[j] = *max_element(vCstr_j.begin(), vCstr_j.end());
	}


	return maxCstrViol;
}

int selectionRule(double xi_fit, vector<double> xi_cstrViol, double ui_fit, vector<double> ui_cstrViol, vector<double> cstrViolMax_G);

vector<double> binomialCrossover(const vector<double>& xi, const vector<double>& vi, const double& CR);
vector<double> exponentialCrossover(const vector<double>& xi, const vector<double>& vi, const double& CR);
vector<double> arithmeticCrossover_line(const vector<double>& xi, const vector<double>& vi, const double& CR);
vector<double> arithmeticCrossover_square(const vector<double>& xi, const vector<double>& vi, const double& CR);

void enforceBoxCstr(vector<double>& x, const vector<double>& lb, const vector<double>& ub);
void enforceBoxCstr(vector<double>& x, const vector<double>& lb, const vector<double>& ub, const vector<int>& boxCstrMode);

template<typename T_agent>
T_agent createMutantAgent(const vector<T_agent>& pop, const T_agent& Gbest, const int& iAgent, const T_agent a_xi, const double& F, const int& iStrategy = 2)
{
	/// T_agent � un template che assume T = C_agent
	///  ma va bene anche per C_jagent
	///
	/// operate on genotipe only:
	/// fit, cstrViol:	not evaluated
	/// F, CR, iStrategy: same as  "a_xi"
	///


	T_agent a_vi = a_xi;

	int np = pop.size();
	int d = a_xi.x.size();



	if (iStrategy == 1)
	{
		/// DE_rand_1:  (3)
		///		vi = xr1 + F*(xr2 - xr3);

		int nRandAgent = 3;  //n. of random agents required for this mutation strategy		
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
			a_vi.x[j] = pop[irr[0]].x[j] + F*(pop[irr[1]].x[j] - pop[irr[2]].x[j]);
	}
	else if (iStrategy == 2)
	{
		/// DE_best_1:  (2)
		///		vi = x_Gbest + F*(xr1-xr2);

		int nRandAgent = 2;  //n. of random agents required for this mutation strategy
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
			a_vi.x[j] = Gbest.x[j] + F*(pop[irr[0]].x[j] - pop[irr[1]].x[j]);
	}
	else if (iStrategy == 3)
	{
		/// DE_current2best_1: (2)
		///		vi = xi + F*(x_Gbest - xi) + F*(xr1 - xr2);

		int nRandAgent = 2;  //n. of random agents required for this mutation strategy
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
			a_vi.x[j] = a_xi.x[j] + F*(Gbest.x[j] - a_xi.x[j]) + F*(pop[irr[0]].x[j] - pop[irr[1]].x[j]);

	}
	else if (iStrategy == 4)
	{

		/// DE_best_2:   (4)
		///		vi = x_Gbest + F*(xr1 - xr2) + F*(xr3 - xr4);

		int nRandAgent = 4;  //n. of random agents required for this mutation strategy
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
		{
			a_vi.x[j] = Gbest.x[j] + F*(pop[irr[0]].x[j] - pop[irr[1]].x[j]) + F*(pop[irr[2]].x[j] - pop[irr[3]].x[j]);
		}

	}
	else if (iStrategy == 5)
	{
		/// DE_rand_2:   (5)
		///		vi = xr1 + F*(xr2 - xr3) + F*(xr4 - xr5);

		int nRandAgent = 5;  //n. of random agents required for this mutation strategy
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
		{
			a_vi.x[j] = pop[irr[0]].x[j] + F*(pop[irr[1]].x[j] - pop[irr[2]].x[j]) + F*(pop[irr[3]].x[j] - pop[irr[4]].x[j]);
		}
	}
	else if (iStrategy == 6)
	{
		/// GA_like_1:   (1)
		///		vi = xi + alfa*( xi - xr1) 		alfa ~ U(0, 1.10)
		const double alfaMax = 1.1;
		const double alfaMin = 0.;

		int nRandAgent = 1;  //n. of random agents required for this mutation strategy
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
		{
			double alfa = alfaMin + (double)(mygenrand()*(alfaMax - alfaMin));
			a_vi.x[j] = a_xi.x[j] + alfa*(a_xi.x[j] - pop[irr[0]].x[j]);
		}
	}
	else if (iStrategy == 7)
	{
		/// GA_like_2:   (2)
		///		vi = xr1 + alfa*( xr2 - xr1) 		alfa ~ U(0, 1.10)
		const double alfaMax = 1.1;
		const double alfaMin = 0.;

		int nRandAgent = 2;  //n. of random agents required for this mutation strategy
		vector<int> irr = i_rrSelect(nRandAgent, np, iAgent);

		for (int j = 0; j < d; j++)
		{
			double alfa = alfaMin + (double)(mygenrand()*(alfaMax - alfaMin));
			a_vi.x[j] = pop[irr[0]].x[j] + alfa*(pop[irr[1]].x[j] - pop[irr[0]].x[j]);
		}
	}


	return a_vi;
}

void rndAdaptRule_F(double& F);
void rndAdaptRule_CR(double& CR);
void rndAdaptRule_S(int& iStrategy);
void ditherAdaptRule_S(int& iStrategy);
void uniformRouletteRule_S(int& iStrategy);	//double freqDEmutationMode[] = { 2., 1., 1., 1., 2. };
void selfAdaptRule(double& F, double& CR, int& iStrategy);

void saveOutputOnDisk(ostream& os, const int& nGen, const C_agent& Gbest);
void saveOutputOnDisk(ostream& os, const int& nGen, const C_agent& Gbest, const double& fitMean, const double& fitStdev);

void evalPopFitStat(const vector<C_jagent>& pop, double& fitMean, double& fitStdev);

void pruneByClustering(const double& rho, vector<C_jagent>& pop, vector<double>& lb, vector<double>& ub, vector<int> boxCstrMode);


vector<C_jagent> popReduction(vector<C_jagent> pop, int& np);


void dispInit(ostream& os);
void dispAdapt(ostream& os, const int& nGen, const vector<C_jagent>& pop);




double epsRule(const double& eps0, const int& t, const int& Tc);

//// Efficiently sampling points uniformly from the surface of an n-sphere
//#include <random>
//void  uniSphere(int n, vector<double>& y)
//{
//
//	std::random_device rd;
//	std::mt19937 gen(rd());
//
//	std::normal_distribution<> normalDistr(0, 1.);
//	// Y is the output direction
//	vector<double> x(n);
//	for (int i = 0; i < n; i++)
//		x[i] = normalDistr(gen);
//
//	double sumSquared = 0;
//	for (int i = 0; i < n; i++)
//		sumSquared += x[i] * x[i];
//	sumSquared = sqrt(sumSquared);
//
//	for (int i = 0; i < n; i++)
//		y[i] = x[i] / sumSquared;
//
//};




 