/*!

versione 1.2
requires: "c_agent.h" v1.0 +

++ modificato  - MinCstrProblemAbstract:
aggiunto "vector<int> boxCstrMode" to allow circular box constraints support




*/
#include <assert.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


#include "c_agent.h"


using namespace std;

#ifndef JEVAL_H
#define JEVAL_H


const double hCstrTol = 0.0001; //Tolleranza sui vincoli di uguaglianza

extern int iTest;



extern int FES; //n� of Fitness evaluation


#endif // !JEVAL_H



#ifndef EvaluateOp_h
#define EvaluateOp_h


class EvaluateOp
{
public:
	virtual double evaluateFitness(const vector<double>& x) const { return -1; };
	virtual vector<double> evaluateCstr(const vector<double>& x) const { return vector<double>(1, -1); };

	virtual void eval_Fit_and_Cstr(C_agent& agent) const {
		agent.fit = evaluateFitness(agent.x);
		agent.cstrViol = evaluateCstr(agent.x);
	};
};
typedef EvaluateOp* EvaluateOpP;	//puntatore alla classe EvaluateOp


class MinCstrProblemAbstract : public EvaluateOp
{
public:
	int d;
	int nCstr;
	vector<double> lb, ub;
	vector<int> boxCstrMode;

	MinCstrProblemAbstract(int d_in=1, int nCstr_in=1) : \
		d(d_in), nCstr(nCstr_in) {
			lb = vector<double>(d, -1e6);
			ub = vector<double>(d, +1e6);
			boxCstrMode = vector<int>(d, 0);	//default = classical box
		};

	MinCstrProblemAbstract(int nCstr_in, vector<double> lb_in, vector<double> ub_in) : \
		nCstr(nCstr_in), lb(lb_in), ub(ub_in), d(lb_in.size()) {
			boxCstrMode = vector<int>(d, 0);	//default = classical box
		};

	MinCstrProblemAbstract(int nCstr_in, vector<double> lb_in, vector<double> ub_in, vector<int> boxCstrMode_in) : \
		nCstr(nCstr_in), lb(lb_in), ub(ub_in), d(lb_in.size()), boxCstrMode(boxCstrMode_in) { };


	//virtual functions
	virtual void postproc(vector<double>& xOpt, ostream& os) const {
		cout << "default postprocess: do nothing" << endl;
		//cin.get();
		return;
	};
};


typedef  MinCstrProblemAbstract* MinCstrPblmbAbsP;





// class MyProb : public MinCstrProblemAbstract
// {
// 	/*
// 	Esempio:
// 	risolvi il problema parametrico
// 	min(x) : (x-p)^2 + y^2
// 	-5<x<5
// 	*/

// private:
// 	static const int _d = 2;
// 	static const int _nCstr = 1;


// public:

// 	double p;	//parametro del problema
// 	//
// 	double v1;

// 	MyProb() : MinCstrProblemAbstract(_d, _nCstr), p(0) {};
// 	MyProb(double parametro) : MinCstrProblemAbstract(_d, _nCstr),	p(parametro) 
// 	{
// 		lb = vector<double>{-5., -5.};
// 		ub = vector<double>{5., 5.};
// 	};

// 	double evaluateFitness(const vector<double>& x) const {
// 		return x[0];
// 		// cout << "x = " << x[0] << '\t' << x[1] << endl;
// 		// cin.get();

// 	};
// 	vector<double> evaluateCstr(const vector<double>& x) const {
// 		double g[_nCstr];
// 		g[0] = pow(x[0] - p, 2) + pow(x[1], 2);

// 		vector<double> 	cstrViol;
// 		cstrViol.push_back(g[0]);

// 		return cstrViol;
// 	};

// };



inline double inequalityConstraint(const double& g, const double& eps);
inline double equalityConstraint(const double& h, const double& eps);

#endif // EvaluateOp_h