#include <iostream>
#include "FileDealing.hpp"
#include "HeisenbergCanvas.hpp"

const int BIN_NUM = 30;
const int SWEEPS_PER_BIN = 100;
const int THERMALIZE_SWEEPS = 1000;
std::string FILE_NAME = "BinderData.txt";
std::string HEADER = "BINDER RATIO";




int main()
{
    int L = 10;
    double T = 5;
    std::cout << "The temperature is " << T << std::endl;
    HeisenbergCanvas h(L,T);
    h.Initialize();

    std::ofstream fout;
    WriteHeader(fout, FILE_NAME, HEADER);
    
    // Thermalize, Metro : Wolff = 9 : 1
    for (int walk=0; walk<THERMALIZE_SWEEPS; walk++)
    {
        for (int i=0;i<9;++i) h.MetropolisWalk();
        h.LBWolff();
    }
    
    std::cout << h.CalculateM() << std::endl;

    // Observe
    std::vector<double> susceptibility_set;
    for (int bin=0;bin<BIN_NUM;++bin)
    {
        h.CountingSuscept(SWEEPS_PER_BIN, susceptibility_set);
    }
    
    h.OutStream(fout, FILE_NAME, susceptibility_set);
    
    std::cout << "Simulation Success!" << std::endl;
    system("pause");

}

