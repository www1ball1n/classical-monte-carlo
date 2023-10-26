#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <ctime>

const int L = 40;
const int J = 1;
const int BIN_NUM = 30;
const int SWEEPS_PER_BIN = 500;
const int THERMALIZE_SWEEPS = 4000;

using Matrix = std::array<std::array<int,L>,L>;


double UniformRandom()
{
    return rand() / (RAND_MAX + 1.); 
}


void Initialize(Matrix& canvas)
{
    for (size_t i = 0; i < L; i++)
        for (size_t j = 0; j < L; j++)
                canvas[i][j] = UniformRandom() < 0.5 ? 1 : -1;
    return;
}


int GetMagnetization(const Matrix& canvas)
{
    int sum = 0;
    for (size_t i = 0; i < L; i++)
        for (size_t j = 0; j < L; j++)
            sum += canvas[i][j];
    return sum;
}

int P_Plus(int i)
{
    return (i + 1) % L;
}

int P_Minus(int i)
{
    return (i + L - 1) % L;
}

double LocalEnergy(const Matrix& canvas, int i, int j)
{
    int spin = canvas[i][j];
    int surroundings = canvas[P_Plus(i)][j] + canvas[P_Minus(i)][j]
                     + canvas[i][P_Plus(j)] + canvas[i][P_Minus(j)];
    return -1.0 * spin * surroundings;
                
}

void MetropolisFlip(Matrix& canvas, double T, int i, int j)
{
    double energy_difference = -2 * LocalEnergy(canvas, i, j);
    double probability = std::exp(-energy_difference / T);
    if (energy_difference < 0 || UniformRandom() < probability)
    {
        canvas[i][j] = -1 * canvas[i][j];
        return;
    }
    else
        return;
}


void MetropolisSweep(Matrix& canvas, double T)
{
    for (size_t i = 0; i < L; i++)
        for (size_t j = 0; j < L; j++)
            MetropolisFlip(canvas, T, i, j);
}


void SweepLoop(Matrix& canvas, int number, double T)
{
    for (size_t i = 0; i < number; i++)
        MetropolisSweep(canvas, T);
}

void CountingCapacity(Matrix& canvas, int number, double T, std::vector<double>& Xset)
{
    long long m_sum = 0;
    long long m2_sum = 0;
    int m = 0;
    for (size_t i = 0; i < number; i++)
    {
        MetropolisSweep(canvas, T);
        m = GetMagnetization(canvas);
        m_sum += m;
        m2_sum += m * m;
    }
    m_sum /= number;
    m2_sum /= number;
    double susceptibility = (m2_sum - m_sum * m_sum) / (T * std::pow(L,4));
    Xset.push_back(susceptibility);    
}

int main()
{
    Matrix canvas;
    srand((unsigned)time(NULL));

    std::ofstream fout;
    fout.open("data.txt", std::ios::out);
    fout << "2023/10/25 ISING SUSCEPTIBILITY" << std::endl; 
    fout.close();

    for (double T = 1.55; T < 3.45; T += 0.05)
    {
        std::cout << "The temperature is " << T << std::endl;
        Initialize(canvas);
        SweepLoop(canvas, THERMALIZE_SWEEPS, T); // Thermalize

        std::vector<double> susceptibility_set;
        for (int bin = 0; bin < BIN_NUM; bin++)
            CountingCapacity(canvas, SWEEPS_PER_BIN, T, susceptibility_set);

        // file outstream
        fout.open("data.txt", std::ios::out|std::ios::app);
        fout << T << std::endl; 
        for (int i = 0; i < susceptibility_set.size(); i++)
            fout << susceptibility_set[i] << std::endl;
        fout.close();

        std::cout << "The Magnetization is " << 1.0 * std::abs(GetMagnetization(canvas)) / L / L << std::endl;
    }
    
    system("pause");
}

