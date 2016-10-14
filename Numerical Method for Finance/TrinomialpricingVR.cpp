/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TrinomialpricingVR.h"
#include <vector>
#include <iostream>
#include "BlackScholes.h"
#include <math.h>

TrinoTreeVR::TrinoTreeVR(double m, double s, double k,
        double sig, double Q, double R,
        bool call, bool euro) : maturity(m), S(s), K(k), sigma(sig), q(Q), r(R), isCall(call), isEuro(euro) {

}

void TrinoTreeVR::pickType(bool e, bool c) {
    isEuro = e;
    isCall = c;
}

double TrinoTreeVR::SolveTT(int n) {
    double deltaT = maturity / n;
    double up = exp(sigma * sqrt(3 * deltaT));
    double pu = (1. / 6. + (r - q - pow(sigma, 2) / 2) * sqrt(deltaT / 12 / pow(sigma, 2))) * exp(-r * deltaT);
    double pm = 2. / 3. * exp(-r * deltaT);
    double pd = exp(-r * deltaT) - pu - pm;
    BlackScholes BS(S, K, maturity, sigma, r, q);
    std::vector<double> price(2 * n + 1, 0);
    std::vector<double> euroTree;

    for (int i = 0; i <= 2 * n; i++) {
        if (!isCall)
            price[i] = K - S * pow(up, n - i);
        else
            price[i] = S * pow(up, n - i) - K;
        if (price[i] < 0)
            price[i] = 0;
        if (!isCall && !isEuro)
            euroTree.push_back(price[i]);
    }

    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j <= 2 * i; j++) {
            price[j] = pu * price[j] + pm * price[j + 1] + pd * price[j + 2];
            if (!isCall && !isEuro) {
                euroTree[j] = pu * euroTree[j] + pm * euroTree[j + 1] + pd * euroTree[j + 2];
                double exercise = K - S * pow(up, i - j);
                if (price[j] < exercise)
                    price[j] = exercise;
//                BlackScholes temp(S * pow(up, i - j), K, deltaT * (n - i), sigma, r, q);
//                price[j] = price[j] + temp.PutPrice() - euroTree[j];
            }
        }
        if (i == 1){
            delta = (price[0] - price[2]) / S / (up - 1 / up) - 0.378673683
                    -(euroTree[0] - euroTree[2]) / S / (up - 1 / up);
        }

        if (i == 2){
            gamma = ((price[0] - price[2]) / S / (up * up - 1)-(price[2] - price[4]) / S / (1 - 1 / up / up)) / S / (up - 1 / up)
                    -((euroTree[0] - euroTree[2]) / S / (up * up - 1)-(euroTree[2] - euroTree[4]) / S / (1 - 1 / up / up)) / S / (up - 1 / up)
                    + 0.030708035;
        }
    }
    theta = (price[1] - price[0]) / deltaT - 1.895422173 - (euroTree[1] - euroTree[0]) / deltaT;
    if(!isCall && !isEuro)
        price[0] = price[0] + BlackScholes(S, K, maturity, sigma, r, q).PutPrice() - euroTree[0];
//    std::cout << (price[1] - price[0]) / deltaT;
//    std::cout << price[0] << " " << price[1] << " " << price[2] << " ";
    return price[0];
}

double TrinoTreeVR::SolveATT(int n) {
    return (SolveTT(n) + SolveTT(n - 1)) / 2;
}

double TrinoTreeVR::SolveBTT(int n) {
    double deltaT = maturity / n;
    double up = exp(sigma * sqrt(3 * deltaT));
    double pu = (1. / 6. + (r - q - pow(sigma, 2) / 2) * sqrt(deltaT / 12 / pow(sigma, 2))) * exp(-r * deltaT);
    double pm = 2. / 3. * exp(-r * deltaT);
    double pd = exp(-r * deltaT) - pu - pm;
    std::vector<double> price(2 * n + 1, 0);
    std::vector<double> euroTree;

    for (int i = 0; i <= (2 * n - 2); i++) {
        BlackScholes temp(S * pow(up, n - 1 - i), K, deltaT, sigma, r, q);
        if (!isCall)
            price[i] = std::max(temp.PutPrice(), K - S * pow(up, n - 1 - i));
        else
            price[i] = temp.CallPrice();

        if (!isEuro && !isCall)
            euroTree.push_back(BlackScholes(S * pow(up, n - 1 - i), K, deltaT, sigma, r, q).PutPrice());
    }

    for (int i = n - 2; i >= 0; i--) {
        for (int j = 0; j <= 2 * i; j++) {
            price[j] = (pu * price[j] + pm * price[j + 1] + pd * price[j + 2]);
            if (!isCall && !isEuro) {
                euroTree[j] = (pu * euroTree[j] + pm * euroTree[j + 1] + pd * euroTree[j + 2]);
                double exercise = K - S * pow(up, i - j);
//                BlackScholes temp(S * pow(up, i - j), K, deltaT * (n - i), sigma, r, q);
                if (price[j] < exercise)
                    price[j] = exercise;
//                price[j] = price[j] + temp.PutPrice() - euroTree[j];
            }
        }
        if (i == 1)
            delta = (price[0] - price[2]) / S / (up - 1 / up) - 0.378673683
                    -(euroTree[0] - euroTree[2]) / S / (up - 1 / up);

        if (i == 2){
            gamma = ((price[0] - price[2]) / S / (up * up - 1)-(price[2] - price[4]) / S / (1 - 1 / up / up)) / S / (up - 1 / up)
                    -((euroTree[0] - euroTree[2]) / S / (up * up - 1)-(euroTree[2] - euroTree[4]) / S / (1 - 1 / up / up)) / S / (up - 1 / up)
                    + 0.030708035;
//            gamma = ((euroTree[0] - euroTree[2]) / S / (up * up - 1)-(euroTree[2] - euroTree[4]) / S / (1 - 1 / up / up)) / S / (up - 1 / up);
        }
    }
//    std::cout << gamma;
    theta = (price[1] - price[0]) / deltaT - 1.895422173 - (euroTree[1] - euroTree[0]) / deltaT;
    if(!isCall && !isEuro)
        price[0] = price[0] + BlackScholes(S, K, maturity, sigma, r, q).PutPrice() - euroTree[0];
    return price[0];
}

double TrinoTreeVR::SolveBTTR(int n) {
    double temp1 = SolveBTT(n);
    double temp2 = delta;
    double temp3 = gamma;
    double temp4 = theta;
    double temp5 = SolveBTT(n / 2);
    delta = 2 * temp2 - delta;
    gamma = 2 * temp3 - gamma;
    theta = 2 * temp4 - theta;
    return 2 * temp1 - temp5;
}

double TrinoTreeVR::Delta(int n, int i) {
    switch (i) {
        case 1: SolveTT(n);
            break;
        case 2: SolveATT(n);
            break;
        case 3: SolveBTT(n);
            break;
        case 4: SolveBTTR(n);
            break;
    }
    return delta;
}

double TrinoTreeVR::Gamma(int n, int i) {
    switch (i) {
        case 1: SolveTT(n);
            break;
        case 2: SolveATT(n);
            break;
        case 3: SolveBTT(n);
            break;
        case 4: SolveBTTR(n);
            break;
    }
    return gamma;
}

double TrinoTreeVR::Theta(int n, int i) {
    switch (i) {
        case 1: SolveTT(n);
            break;
        case 2: SolveATT(n);
            break;
        case 3: SolveBTT(n);
            break;
        case 4: SolveBTTR(n);
            break;
    }
    return theta;
}

TrinoTreeVR::TrinoTreeVR(const TrinoTreeVR& orig) {
}

TrinoTreeVR::~TrinoTreeVR() {
}