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
    void CountingSuscept(int number, std::vector<double>& Xset);
    void OutStream(std::ofstream& fout,std::string filename, std::vector<double> & sset);
    void BalanceTest(int number);
    void CountingM(int number, std::vector<double>& Xset);
    void MixWalk();
    
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
    
    do
    {
        // 获取x,y坐标
        currentSpin = stack.back();
        x = currentSpin[0];
        y = currentSpin[1];
        stack.pop_back();

        canvas[x][y] = canvas[x][y] - 2 * mirror * Dot(mirror,canvas[x][y]); // Flip
        canvas[x][y] = canvas[x][y] * (1.0/norm(canvas[x][y]));
        isingCanvas[x][y] *= -1;
        
        double cumulative[L/2-1+1];
        double sumTemp;

        // Right
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2-1;i++)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[(L+x+i)%L][y] / (i * i * T)); 
            cumulative[i] = 1 - std::exp(sumTemp); 
        }

        int k = 0;
        // 每次的conditionalCumulative都是chainLength+1的大小
        while (true)
        {
            double dice = rd.UniformRandom();
            dice = cumulative[k]+(1-cumulative[k]) * dice;
            if (dice >= cumulative[L/2-1]) break; // 构造结束

            k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
            for (int i=1;i<k;i++) isConsidered[(L+x+i)%L][y] = 1;
            int chosenX = (L+x+k)%L, chosenY = y;
            if (!isConsidered[chosenX][chosenY])
            {
                std::array<int,2> coordinate = {chosenX,chosenY};
                stack.push_back(coordinate);
            }
            isConsidered[chosenX][chosenY] = true;
        }

        // Left
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2-1;i++)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[(L+x-i)%L][y] / (i * i * T)); 
            cumulative[i] = 1 - std::exp(sumTemp); 
        }

        k = 0;
        // 每次的conditionalCumulative都是chainLength+1的大小
        while (true)
        {
            double dice = rd.UniformRandom();
            dice = cumulative[k]+(1-cumulative[k]) * dice;
            if (dice >= cumulative[L/2-1]) break; // 构造结束

            k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
            for (int i=1;i<k;i++) isConsidered[(L+x-i)%L][y] = 1;
            int chosenX = (L+x-k)%L, chosenY = y;
            if (!isConsidered[chosenX][chosenY])
            {
                std::array<int,2> coordinate = {chosenX,chosenY};
                stack.push_back(coordinate);
            }
            isConsidered[chosenX][chosenY] = true;
        }

        // Up
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2-1;i++)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[x][(L+y+i)%L] / (i * i * T)); 
            cumulative[i] = 1 - std::exp(sumTemp); 
        }

        k = 0;
        while (true)
        {
            double dice = rd.UniformRandom();
            dice = cumulative[k]+(1-cumulative[k]) * dice;
            if (dice >= cumulative[L/2-1]) break; // 构造结束

            k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
            for (int i=1;i<k;i++) isConsidered[x][(L+y+i)%L] = 1;
            int chosenX = x, chosenY = (L+y+k)%L;
            if (!isConsidered[chosenX][chosenY])
            {
                std::array<int,2> coordinate = {chosenX,chosenY};
                stack.push_back(coordinate);
            }
            isConsidered[chosenX][chosenY] = true;
        }

        // Down
        sumTemp = 0;
        cumulative[0] = 0;
        for (int i=1;i<=L/2-1;i++)
        {
            // 按道理可以用指针运算化简，先放在这里
            sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[x][(L+y-i)%L] / (i * i * T)); 
            cumulative[i] = 1 - std::exp(sumTemp); 
        }

        k = 0;
        while (true)
        {
            double dice = rd.UniformRandom();
            dice = cumulative[k]+(1-cumulative[k]) * dice;
            if (dice >= cumulative[L/2-1]) break; // 构造结束

            k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
            for (int i=1;i<k;i++) isConsidered[x][(L+y-i)%L] = 1;
            int chosenX = x, chosenY = (L+y-k)%L;
            if (!isConsidered[chosenX][chosenY])
            {
                std::array<int,2> coordinate = {chosenX,chosenY};
                stack.push_back(coordinate);
            }
            isConsidered[chosenX][chosenY] = true;
        }

    } while (!stack.empty());

}

void HeisenbergCanvas::LBWalk()
{
    for (int walk=0;walk<10;walk++) LBWolff();
}

void HeisenbergCanvas::MixWalk()
{
    LBWolff();
    for (int i=0;i<5;i++) MetropolisWalk();
}

void HeisenbergCanvas::CountingSuscept(int number, std::vector<double>& Xset)
{
    long double m4_sum = 0;
    long double m2_sum = 0;
    double m = 0;
    for (size_t i = 0; i < number/10; i++)
    {
        for (int j=0;j<10;j++)
        {
            MixWalk();
            m = CalculateM();
            m2_sum += m * m;
            m4_sum += pow(m,4);
        }
    }
    m4_sum /= number;
    m2_sum /= number;
    double susceptibility = m4_sum / (m2_sum * m2_sum);
    Xset.push_back(susceptibility);
    // test
    std::cout << susceptibility << std::endl;
}

void HeisenbergCanvas::CountingM(int number, std::vector<double>& Xset)
{
    long double m_sum = 0;
    double m = 0;
    for (size_t i = 0; i < number/10; i++)
    {
        for (int j=0;j<10;j++)
        {
            MixWalk();
            m = CalculateM();
            m_sum += m;
        }
    }
    m_sum /= number;
    double susceptibility = m_sum;
    Xset.push_back(susceptibility);
    // test
    std::cout << susceptibility << std::endl;
}

void HeisenbergCanvas::BalanceTest(int number)
{
    long double m_sum = 0;
    long double m2_sum = 0;
    double mx;
    for (size_t i = 0; i < number/10; i++)
    {
        for (int j=0;j<9;j++)
        {
            MetropolisWalk();
            mx = canvas[0][0][0];
            m_sum += mx;
            m2_sum += mx * mx;
            
        }
        LBWolff();
        mx = canvas[0][0][0];
        m_sum += mx;
        m2_sum += mx * mx;
    }
    m_sum /= number;
    m2_sum /= number;
    // test
    std::cout << "mx average = " << m_sum << std::endl;
    std::cout << "mx2 average = " << m2_sum << std::endl;
}

