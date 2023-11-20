#pragma once
#include "LuijtenBloete.hpp"
#include "Wolff.hpp"

// 10 times LB
void HeisenbergCanvas::LBWalk()
{
    for (int walk=0;walk<10;walk++) LBWolff();
}

// one LB, then Metropolis sweep L^2
void HeisenbergCanvas::MixWalk()
{
    LBWolff();
    MetropolisWalk();
}

// one Wolff, then Metropolis sweep L^2
void HeisenbergCanvas::MixWolff()
{
    Wolff();
    MetropolisWalk();
}

void HeisenbergCanvas::CountingBinder(int number, std::vector<double>& Xset)
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