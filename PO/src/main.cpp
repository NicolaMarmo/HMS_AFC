#include "../include/C_prb_RR_HMS.h"  
#include "/home/spectral01/MARMO/Tools/tools.h" 

#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <thread>

using namespace std;

void printASCII(string filename){
    string line;
    ifstream infile;
    infile.open(filename);
    if(infile.is_open()){
        while(getline(infile, line)){
            cout << line << endl;
        }
        cout << endl;
    }
    else{
        cout << "File not found" << endl;
    }
    infile.close();
}

int main(){

    printASCII("../JAXA.txt");

    // Hybrid Multi Shooting
    test_RR_HMS();

    // VectorXd s0(6); s0 << 1.0277926091, 0.0, -0.1858044184, 0.0, -0.1154896637, 0.0;
    // C_simulazione sim(6, 2000, 1, 1e-8, 1e-10);
    // sim.SetY0_Sim(s0);
    // EoM_CR3BP EoM(0.012150584269940, 1.5872714606);
    // sim.Start_Sim(EoM);

    // system("cd ..; rm XX.dat");
    // ofstream out_sol("../XX.dat");
    // for (int i=0; i<sim.Y.cols(); i++)
    // {
    //     out_sol << fixed <<setw(12) << setprecision(6) << 1.5872714606*sim.T(i);
    //     for (int j=0; j<sim.Y.rows(); j++)
    //     {
    //         out_sol << fixed << setw(12) << setprecision(6) << sim.Y(j,i);
    //     }
    //     out_sol << endl;
    // }
    // out_sol.close();

    // VectorXd XXf(6); XXf = sim.Y.rightCols(1);
    // VectorXd dX(6); dX = XXf - s0;
    // cout << "dX = " << dX.transpose() << endl;

    // return 0;

}



