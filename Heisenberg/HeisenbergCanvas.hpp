#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "RandomNumber.hpp"
#include "Utils.hpp"

using Spin = std::array<double,3>;

class HeisenbergCanvas
{

public:
    int L;
    double T;
    std::vector<std::vector<Spin>> canvas;
    RandomNumber rd;
    std::vector<double> interaction;
    
public:
    HeisenbergCanvas(int size, double temperature)
    {
        L = size;
        T = temperature;
        canvas = std::vector<std::vector<Spin>>(size, std::vector<Spin>(size));
        rd = RandomNumber(size);
    }

    void Initialize();
    double CalculateM();
    Spin SurroundingOf(int x, int y);
    void MetropolisWalk();
    void LBWolff();
    void LBWalk();
    void CountingBinder(int number, std::vector<double>& Xset);
    void OutStream(std::ofstream& fout,std::string filename, std::vector<double> & sset);
    void BalanceTest(int number);
    void CountingM(int number, std::vector<double>& Xset);
    void MixWalk();
    void Wolff();
    void MixWolff();
    
};

void HeisenbergCanvas::Initialize()
{
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            rd.RandomizeSpin(canvas[i][j]);
        }
    }
    
    interaction.push_back(0.0);
    for (int i=1;i<=L/2;i++)
    {
        double sum = pow(M_PI / (L * sin(i * M_PI / L)),2);
        // for (int j=0;j<10;j++) sum += 1.0 / pow(i + j * L, 2);
        interaction.push_back(sum / T);
    }

    for (int i=0;i<=L/2;i++)
    {
        std::cout << interaction[i] << " ";
    }
    std::cout << std::endl;
}

double HeisenbergCanvas::CalculateM()
{
    Spin sum = {0,0,0};
    int count[2] = {1,-1};
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            sum += canvas[i][j] * count[(i+j)%2];
        }
    }
    return (norm(sum) / (L * L));
}

Spin HeisenbergCanvas::SurroundingOf(int x, int y)
{
    Spin res = {0,0,0};
    for (int i=1;i<L/2;i++)
    {
        res += canvas[(L+x+i)%L][y] * interaction[i];
        res += canvas[(L+x-i)%L][y] * interaction[i];
        res += canvas[x][(L+y-i)%L] * interaction[i];
        res += canvas[x][(L+y-i)%L] * interaction[i];
    }
    // We ensure L is even
    res += canvas[(L/2+x)%L][y] * (4.0/(L*L));
    res += canvas[x][(L/2+x)%L] * (4.0/(L*L));
    return res;
}

void HeisenbergCanvas::MetropolisWalk()
{
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) {
        Spin newSpin = rd.GetRandomSpin();
        double energyDifference = Dot((newSpin - canvas[i][j]),SurroundingOf(i,j));
        if (energyDifference < 0 || rd.UniformRandom() < exp(-energyDifference))
        {
            canvas[i][j] = newSpin;
        }
    }}
}




