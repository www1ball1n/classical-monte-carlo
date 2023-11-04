#pragma once
#include <array>
#include <iostream>
#include <cmath>

using Spin = std::array<double,3>;

// Overload += for an array
Spin& operator+=(Spin &s1, const Spin& s2)
{
    for (int i=0;i<3;i++) s1[i] += s2[i];
    return s1;
}

// operator+ creates a new array.
Spin operator+(const Spin &s1, const Spin& s2)
{
    Spin res;
    for (int i=0;i<3;i++) res[i] = s1[i]+s2[i];
    return res;
}

Spin operator-(const Spin &s1, const Spin& s2)
{
    Spin res;
    for (int i=0;i<3;i++) res[i] = s1[i]-s2[i];
    return res;
}

// operator* creates a new array.
Spin operator* (const Spin& s1, const double c)
{
    Spin res;
    for (int i=0;i<3;i++) res[i] = c * s1[i];
    return res;
} 

Spin operator* (const double c, const Spin& s1)
{
    Spin res;
    for (int i=0;i<3;i++) res[i] = c * s1[i];
    return res;
}

// Calculate norm
double norm(Spin s)
{
    double sum = 0;
    for (int i=0;i<3;i++) sum += s[i] * s[i];
    return sqrt(sum);
}

// cout an array
std::ostream& operator<<(std::ostream& c, const Spin& s)
{
    std::cout << s[0] << " " << s[1] << " " << s[2] << " ";
    return std::cout;
}

// Dot Product
double Dot(const Spin& s1, const Spin& s2)
{
    double res = 0;
    for (int i=0;i<3;i++)
    {
        res += s1[i] * s2[i];
    }
    return res;
}


