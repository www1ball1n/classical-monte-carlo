#pragma once
#include <array>
#include <cmath>
#include <random>
#include <chrono>


class RandomNumber
{
private:
    /* data */
public:
    int L;
    std::default_random_engine generator;
    std::normal_distribution<double> gaussian;
    std::uniform_real_distribution<double> uniform;
    std::uniform_int_distribution<int> uni;
    
    double s0;
    int randint;
    std::array<double,3> s;
public:

    RandomNumber(int size)
    {
        L = size;
        generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
        uni = std::uniform_int_distribution<int>(0,L-1);
    }
    
    // Do not remove default constructor
    RandomNumber()
    {
    
    }

    double& UniformRandom()
    {
        s0 = uniform(generator);
        return s0;
    }

    void RandomizeSpin(std::array<double,3>& spin)
    {
        double sum = 0;
        for (int i=0;i<3;i++) spin[i] = gaussian(generator);
        for (int i=0;i<3;i++) sum += spin[i] * spin[i];
        sum = sqrt(sum);
        for (int i=0;i<3;i++) spin[i] /= sum;
    }

    std::array<double,3>& GetRandomSpin()
    {
        double sum = 0;
        for (int i=0;i<3;i++) s[i] = gaussian(generator);
        for (int i=0;i<3;i++) sum += s[i] * s[i];
        sum = sqrt(sum);
        for (int i=0;i<3;i++) s[i] /= sum;
        return s;
    }

    int& GetRandomInt()
    {
        randint = uni(generator);
        return randint;
    }
};

