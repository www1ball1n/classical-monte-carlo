#pragma once
#include "HeisenbergCanvas.hpp"

void HeisenbergCanvas::LBWolff()
{
    // Spin mirror = rd.GetRandomSpin();
    Spin mirror = {0,0,1};
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
    int clusterSize = 0;
    
    do
    {
        // 获取x,y坐标
        currentSpin = stack.back();
        x = currentSpin[0];
        y = currentSpin[1];
        stack.pop_back();
        clusterSize++;

        Spin spinValue = canvas[x][y];
        canvas[x][y] = spinValue - 2 * mirror * Dot(mirror,spinValue); // Flip
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
            sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[(L+x+i)%L][y] * interaction[i]); 
            cumulative[i] = 1 - std::exp(sumTemp); 
        }

        int k = 0;
        // 每次的Cumulative都是chainLength+1的大小
        while (true)
        {
            double dice = rd.UniformRandom();
            dice = cumulative[k]+(1-cumulative[k]) * dice;
            if (dice >= cumulative[L/2-1]) break; // 构造结束

            k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
            // for (int i=1;i<k;i++) isConsidered[(L+x+i)%L][y] = 1;
            int chosenX = (L+x+k)%L, chosenY = y;
            if (!isConsidered[chosenX][chosenY])
            {
                std::array<int,2> coordinate = {chosenX,chosenY};
                stack.push_back(coordinate);
            }
            isConsidered[chosenX][chosenY] = true;
        }

        // // Left
        // sumTemp = 0;
        // cumulative[0] = 0;
        // for (int i=1;i<=L/2-1;i++)
        // {
        //     // 按道理可以用指针运算化简，先放在这里
        //     sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[(L+x-i)%L][y] / (i * i * T)); 
        //     cumulative[i] = 1 - std::exp(sumTemp); 
        // }

        // k = 0;
        // // 每次的conditionalCumulative都是chainLength+1的大小
        // while (true)
        // {
        //     double dice = rd.UniformRandom();
        //     dice = cumulative[k]+(1-cumulative[k]) * dice;
        //     if (dice >= cumulative[L/2-1]) break; // 构造结束

        //     k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
        //     for (int i=1;i<k;i++) isConsidered[(L+x-i)%L][y] = 1;
        //     int chosenX = (L+x-k)%L, chosenY = y;
        //     if (!isConsidered[chosenX][chosenY])
        //     {
        //         std::array<int,2> coordinate = {chosenX,chosenY};
        //         stack.push_back(coordinate);
        //     }
        //     isConsidered[chosenX][chosenY] = true;
        // }

        // Up
        // sumTemp = 0;
        // cumulative[0] = 0;
        // for (int i=1;i<=L/2-1;i++)
        // {
        //     // 按道理可以用指针运算化简，先放在这里
        //     sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[x][(L+y+i)%L] * interaction[i]); 
        //     cumulative[i] = 1 - std::exp(sumTemp); 
        // }

        // k = 0;
        // while (true)
        // {
        //     double dice = rd.UniformRandom();
        //     dice = cumulative[k]+(1-cumulative[k]) * dice;
        //     if (dice >= cumulative[L/2-1]) break; // 构造结束

        //     k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
        //     for (int i=1;i<k;i++) isConsidered[x][(L+y+i)%L] = 1;
        //     int chosenX = x, chosenY = (L+y+k)%L;
        //     if (!isConsidered[chosenX][chosenY])
        //     {
        //         std::array<int,2> coordinate = {chosenX,chosenY};
        //         stack.push_back(coordinate);
        //     }
        //     isConsidered[chosenX][chosenY] = true;
        // }

        // // Down
        // sumTemp = 0;
        // cumulative[0] = 0;
        // for (int i=1;i<=L/2-1;i++)
        // {
        //     // 按道理可以用指针运算化简，先放在这里
        //     sumTemp += std::min(0.0, -2 * isingCanvas[x][y] * isingCanvas[x][(L+y-i)%L] / (i * i * T)); 
        //     cumulative[i] = 1 - std::exp(sumTemp); 
        // }

        // k = 0;
        // while (true)
        // {
        //     double dice = rd.UniformRandom();
        //     dice = cumulative[k]+(1-cumulative[k]) * dice;
        //     if (dice >= cumulative[L/2-1]) break; // 构造结束

        //     k = std::upper_bound(cumulative,cumulative+L/2-1+1,dice)-cumulative;// 查找首个
        //     for (int i=1;i<k;i++) isConsidered[x][(L+y-i)%L] = 1;
        //     int chosenX = x, chosenY = (L+y-k)%L;
        //     if (!isConsidered[chosenX][chosenY])
        //     {
        //         std::array<int,2> coordinate = {chosenX,chosenY};
        //         stack.push_back(coordinate);
        //     }
        //     isConsidered[chosenX][chosenY] = true;
        // }

    } while (!stack.empty());

}