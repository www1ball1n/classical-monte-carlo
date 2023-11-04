#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
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
    void CountingSuscept(int number, std::vector<double>& Xset);
    void OutStream(std::ofstream& fout,std::string filename, std::vector<double> & sset);
    
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
}

double HeisenbergCanvas::CalculateM()
{
    Spin sum = {0,0,0};
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            sum += canvas[i][j] * (2*((i+j)%2)-1);
        }
    }
    return (norm(sum) / (L * L));
}

Spin HeisenbergCanvas::SurroundingOf(int x, int y)
{
    Spin res = {0,0,0};
    for (int i=1;i<L/2;i++)
    {
        res += canvas[(L+x+i)%L][y] * (1.0/(i*i));
        res += canvas[(L+x-i)%L][y] * (1.0/(i*i));
        res += canvas[x][(L+y-i)%L] * (1.0/(i*i));
        res += canvas[x][(L+y-i)%L] * (1.0/(i*i));
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
        if (energyDifference < 0 || rd.UniformRandom() < exp(-energyDifference / T))
        {
            canvas[i][j] = newSpin;
        }
    }}
}

void HeisenbergCanvas::LBWolff()
{
    Spin mirror = rd.GetRandomSpin();
    std::vector<std::vector<double>> isingCanvas;
    std::vector<std::vector<bool>> isConsidered;
    std::vector<std::array<int,2>> stack;
    isConsidered = std::vector<std::vector<bool>>(L, std::vector<bool>(L,0));

    // Project to Ising
    isingCanvas = std::vector<std::vector<double>>(L, std::vector<double>(L,0));
    
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) {
        isingCanvas[i][j] = Dot(mirror,canvas[i][j]);
    }}

    int x = rd.GetRandomInt();
    int y = rd.GetRandomInt();
    std::array<int,2> currentSpin = {x,y};
    isConsidered[x][y] = true;
    stack.push_back(currentSpin);
    
    // 参考qy的写法
    do
    {
        currentSpin = stack.back();
        x = currentSpin[0];
        y = currentSpin[1];
        canvas[x][y] = canvas[x][y] - mirror * Dot(mirror,canvas[x][y]); // Flip
        stack.pop_back();
        int cumulative[L/2+1]; 
        double sumTemp;
        int chainLength;
        int conditionalCumulative[L/2+1];


        // 换一下cumulative储存方式 11/4 16:00
        // Right
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2;++i)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += isingCanvas[x][y] * isingCanvas[(L+x+i)%L][y] / (i * i);
            cumulative[i] = 1 - std::exp(2 * sumTemp); 
        }

        chainLength = L/2;
        // 每次的conditionalCumulative都是chainLength+1的大小
        
        while (true)
        {
            double dice = rd.UniformRandom();
            if (dice > cumulative[chainLength]) break; // 构造结束

            int k = std::lower_bound(cumulative,cumulative+chainLength,dice)-cumulative; ++k;// 查找
            for (int i=1;i<=k;i++) isConsidered[(L+x+i)%L][y] = 1;
            if (!isConsidered[(L+x+k)%L][y] && isingCanvas[x][y]*isingCanvas[(L+x+k)%L][y] < 0)
            {
                std::array<int,2> coordinate = {(L+x+k)%L,y};
                stack.push_back(coordinate);
            }
            chainLength -= k; // new chainLength
            
            for (int i=1;i<=chainLength;++i) // calculate new cumulative
                conditionalCumulative[i-1] = (cumulative[k]-cumulative[k+i])/(1-cumulative[k]);

            for (int i=0;i<chainLength;i++) // substitution
                cumulative[i]=conditionalCumulative[i];
        }


        // Left
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2;++i)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += isingCanvas[x][y] * isingCanvas[(L+x-i)%L][y] / (i * i);
            cumulative[i] = 1 - std::exp(2 * sumTemp); 
        }

        chainLength = L/2;
        // 每次的conditionalCumulative都是chainLength+1的大小
        
        while (true)
        {
            double dice = rd.UniformRandom();
            if (dice > cumulative[chainLength]) break; // 构造结束

            int k = std::lower_bound(cumulative,cumulative+chainLength,dice)-cumulative; ++k;// 查找
            for (int i=1;i<=k;i++) isConsidered[(L+x-i)%L][y] = 1;
            if (!isConsidered[(L+x-k)%L][y] && isingCanvas[x][y]*isingCanvas[(L+x-k)%L][y] < 0)
            {
                std::array<int,2> coordinate = {(L+x-k)%L,y};
                stack.push_back(coordinate);
            }
            chainLength -= k; // new chainLength
            
            for (int i=1;i<=chainLength;++i) // calculate new cumulative
                conditionalCumulative[i-1] = (cumulative[k]-cumulative[k+i])/(1-cumulative[k]);

            for (int i=0;i<chainLength;i++) // substitution
                cumulative[i]=conditionalCumulative[i];
        }


        // Up
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2;++i)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += isingCanvas[x][y] * isingCanvas[x][(y+L+i)%L] / (i * i);
            cumulative[i] = 1 - std::exp(2 * sumTemp); 
        }

        chainLength = L/2;
        // 每次的conditionalCumulative都是chainLength+1的大小
        
        while (true)
        {
            double dice = rd.UniformRandom();
            if (dice > cumulative[chainLength]) break; // 构造结束

            int k = std::lower_bound(cumulative,cumulative+chainLength,dice)-cumulative; ++k;// 查找
            for (int i=1;i<=k;i++) isConsidered[x][(y+L+i)%L] = 1;
            if (!isConsidered[x][(y+L+k)%L] && isingCanvas[x][y]*isingCanvas[x][(y+L+k)%L] < 0)
            {
                std::array<int,2> coordinate = {x,(y+L+k)%L};
                stack.push_back(coordinate);
            }
            chainLength -= k; // new chainLength
            
            for (int i=1;i<=chainLength;++i) // calculate new cumulative
                conditionalCumulative[i-1] = (cumulative[k]-cumulative[k+i])/(1-cumulative[k]);

            for (int i=0;i<chainLength;i++) // substitution
                cumulative[i]=conditionalCumulative[i];
        }


        // Down
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2;++i)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += isingCanvas[x][y] * isingCanvas[x][(y+L-i)%L] / (i * i);
            cumulative[i] = 1 - std::exp(2 * sumTemp); 
        }

        chainLength = L/2;
        // 每次的conditionalCumulative都是chainLength+1的大小
        
        while (true)
        {
            double dice = rd.UniformRandom();
            if (dice > cumulative[chainLength]) break; // 构造结束

            int k = std::lower_bound(cumulative,cumulative+chainLength,dice)-cumulative; ++k;// 查找
            for (int i=1;i<=k;i++) isConsidered[x][(y+L-i)%L] = 1;
            if (!isConsidered[x][(y+L-k)%L] && isingCanvas[x][y]*isingCanvas[x][(y+L-k)%L] < 0)
            {
                std::array<int,2> coordinate = {x,(y+L-k)%L};
                stack.push_back(coordinate);
            }
            chainLength -= k; // new chainLength
            
            for (int i=1;i<=chainLength;++i) // calculate new cumulative
                conditionalCumulative[i-1] = (cumulative[k]-cumulative[k+i])/(1-cumulative[k]);

            for (int i=0;i<chainLength;i++) // substitution
                cumulative[i]=conditionalCumulative[i];
        }

    } while (!stack.empty());


}

void HeisenbergCanvas::CountingSuscept(int number, std::vector<double>& Xset)
{
    long double m4_sum = 0;
    long double m2_sum = 0;
    double m = 0;
    for (size_t i = 0; i < number/10; i++)
    {
        for (int j=0;j<9;j++)
        {
            MetropolisWalk();
            m = CalculateM();
            m2_sum += m * m;
            m4_sum += pow(m,4);
        }
        LBWolff();
        m = CalculateM();
        m2_sum += m * m;
        m4_sum += pow(m,4);
    }
    m4_sum /= number;
    m2_sum /= number;
    double susceptibility = m4_sum / (m2_sum * m2_sum);
    Xset.push_back(susceptibility);
    // test
    std::cout << susceptibility << std::endl;
}

