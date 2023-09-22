#include "c_agent.h"
#include <algorithm>
#include "jeval.h"

/***

23/11/2018 - add simple "save to disk" 

*/




#ifndef C_DEME
#define C_DEME




typedef bool *(cmpOpAgents_fptr)(C_agent a_left, C_agent a_right);

class C_deme
{
private:



public:
	int np;
	int nGen = 0;	//deme's age
	int genStall_xGbest = 0;		//n. of stall generation @ x
	int genStall_fitGbest = 0;		//n. of stall generation @ fit
	vector<C_jagent> pop;
	vector<double> maxCstrViol;
	//vector<C_jagent>::iterator it_Gbest;
	C_jagent Gbest;


	//Constructor
	C_deme();

	// crea la popolazione (solo i genotipi)
	C_deme(int np_in, vector<double> lb, vector<double> ub, int creationMethod = 0);



	//memeber functions

	void initializeRndANDSeeded(double soluzioneSeeded, vector<double> lb, vector<double> ub);
	void init_genotype_unifrnd(const vector<double>& lb, const vector<double>& ub);
	//void eval(MinCstrPblmbAbsP ptr_Problem);
	void findBestOfDeme(int _selectionType = 1);


	//virtual void evolve();

	void save(ostream &os) const;




};






#endif // !C_DEME