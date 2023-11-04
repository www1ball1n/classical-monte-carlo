#include <iostream>
#include "HeisenbergCanvas.hpp"


const int BIN_NUM = 30;
const int SWEEPS_PER_BIN = 1000;
const int THERMALIZE_SWEEPS = 5000;



int main()
{
    int L;
    double START_TEMP, END_TEMP, STEP;

    std::cout << "--------------------------------" << "\n";
    std::cout << "DETAILED BALANCE TEST CODE" << "\n";
    std::cout << "--------------------------------" << std::endl;
    
    std::cout << "Please enter the lattice size: ";
    std::cin >> L;
    std::cout << "Please enter the start temperature: ";
    std::cin >> START_TEMP;
    std::cout << "Please enter the end temperature (contained): ";
    std::cin >> END_TEMP;
    std::cout << "Please enter the temperature step: ";
    std::cin >> STEP;


    for (double T=START_TEMP;T<=END_TEMP;T+=STEP)
    {
        std::cout << "The temperature is " << T << std::endl;
        HeisenbergCanvas h(L,T);
        h.Initialize();
    
        // Thermalize, Metro : Wolff = 9 : 1
        for (int walk=0; walk<THERMALIZE_SWEEPS; walk++)
        {
            for (int i=0;i<9;++i) h.MetropolisWalk();
            h.LBWolff();
        }

        // Test
        std::vector<double> susceptibility_set;
        h.BalanceTest(SWEEPS_PER_BIN);
    }
    
    system("pause");
}

