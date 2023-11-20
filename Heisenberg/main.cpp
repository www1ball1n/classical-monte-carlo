#include "Header.hpp"

const int BIN_NUM = 30;
const int SWEEPS_PER_BIN = 2000;
const int THERMALIZE_SWEEPS = 10000;

std::string HEADER = "BINDER RATIO";


int main()
{
    int L;
    std::string FILE_NAME;
    double START_TEMP, END_TEMP, STEP;

    std::cout << "----------------------------------------------------------" << "\n";
    std::cout << "MONTE CARLO SIMULATION FOR SLAC FERMION HEISENBERG MODEL" << "\n";
    std::cout << "----------------------------------------------------------" << std::endl;

    std::cout << "Please enter the lattice size: ";
    std::cin >> L;
    std::cout << "Please enter your file name: ";
    std::cin >> FILE_NAME;
    std::cout << "Please enter the start temperature: ";
    std::cin >> START_TEMP;
    std::cout << "Please enter the end temperature (contained): ";
    std::cin >> END_TEMP;
    std::cout << "Please enter the temperature step: ";
    std::cin >> STEP;

    std::ofstream fout;
    WriteHeader(fout, FILE_NAME, HEADER);

    for (double T=START_TEMP;T<=END_TEMP;T+=STEP)
    {
        std::cout << "The temperature is " << T << std::endl;
        HeisenbergCanvas h(L,T);
        h.Initialize();
    
        for (int walk=0; walk<THERMALIZE_SWEEPS; walk++)
        {
            h.MixWalk();
        }

        // Observe
        std::vector<double> susceptibility_set;
        for (int bin=0;bin<BIN_NUM;++bin)
            h.CountingBinder(SWEEPS_PER_BIN, susceptibility_set); // actually Binder Ratio
        
        h.OutStream(fout, FILE_NAME, susceptibility_set);
    }
    
}

