#include "jDEC_funz_v2.h"
#include "c_deme.h"


#ifndef JDEC_FUNC3_H
#define JDEC_FUNC3_H

#define INNER_OMP_ON 1

//Parametri di Default
const int nTribes_DEFAULT = 4;	//n. di tribù locali
const int nGen_Migration = 100; 
const int nMigrAgent = 1;



//ptr to usrdef postproc function
//typedef PpostProc postproc_fptr



/*   FOR POSTPROCESS PURPOSE   */
const int nGen_ProsProc = 100;
typedef void(*PpostProc)(vector<double>& xGbest, ostream& os);





//forse da spostare
void evalDeme_weightCstrUpdate(C_deme& deme);
void evalDeme(C_deme& deme, MinCstrPblmbAbsP ptr_Problem);
//----------------------------------------------------------//



//void evolveDeme_by_DE(C_deme& deme, MinCstrPblmbAbsP ptr_Problem);
void evolveDeme_by_DE(C_deme& deme, MinCstrPblmbAbsP ptr_Problem, int _adaptModeS = nStrategyMax + 1, int _selectionType = -1);

C_agent multiDeme_optimization(MinCstrPblmbAbsP ptr_Problem, int _np = -1, int _nTribes = -1, int _nGenMAX = -1, int _selectionType = -1);
C_agent multiDeme_optimization_EE(MinCstrPblmbAbsP ptr_Problem, int _np = -1, int _nGenMAX = -1, int _selectionType = -1);


//void initJDEparameters(C_deme& deme);
//void updateJDEparameters(C_deme& deme);
//C_deme createDeme_for_DE(MinCstrPblmbAbsP ptr_Problem, int np);
void initJDEparameters(C_deme& deme, const int& adaptModeF, const int&adaptModeCr, const int& adaptModeS);
void updateJDEparameters(C_deme& deme, const int& adaptModeF, const int&adaptModeCr, const int& adaptModeS);
C_deme createDeme_for_DE(MinCstrPblmbAbsP ptr_Problem, int np, const int& adaptModeF, const int&adaptModeCr, const int& adaptModeS);



void migrationFull(vector<C_deme>::iterator it_tribe0, vector<C_deme>::iterator it_tribe1);
void migration_ring(vector<C_deme>& tribes, int nMigrAgent = 3);
void migration_FW2(vector<C_deme>::iterator it_tribe0, vector<C_deme>::iterator it_tribe1, int nMigrAgent = 3);
void migration_BW2(vector<C_deme>::reverse_iterator it_tribe0, vector<C_deme>::reverse_iterator it_tribe1, int nMigrAgent = 3);


vector<int> linearPartition(const int& npTot, const int& nTribes);
vector<int> powerPartition(const int& npTot, const int& nTribes);





#endif // !JDEC_FUNC3_H


