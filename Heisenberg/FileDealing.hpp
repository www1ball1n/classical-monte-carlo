#pragma once
#include "HeisenbergCanvas.hpp"
#include <fstream>
#include <ctime>
#include <string>

void WriteHeader(std::ofstream& fout, std::string filename, std::string head)
{
    fout.open(filename, std::ios::out);
    fout << head << std::endl; 
    fout.close();
}

void HeisenbergCanvas::OutStream(std::ofstream& fout,std::string filename, std::vector<double> & sset)
{
    fout.open(filename, std::ios::out|std::ios::app);
    fout << T << std::endl; 
    for (int i = 0; i < sset.size(); i++)
        fout << sset[i] << std::endl;
    fout.close();
}