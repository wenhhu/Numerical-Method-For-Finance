/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Trinomialpricing.cpp
 * Author: wenhaohu
 * 
 * Created on October 12, 2016, 4:48 PM
 */

#include "Trinomialpricing.h"
#include <vector>
#include "BlackScholes.h"
#include <math.h>
#include <iostream>

TrinoTree::TrinoTree(double m, double s, double k,
        double sig, double Q, double R,
        bool call, bool euro) : maturity(m), S(s), K(k), sigma(sig), q(Q), r(R), isCall(call), isEuro(euro) {

}

void TrinoTree::pickType(bool e, bool c) {
    isEuro = e;
    isCall = c;
}

double TrinoTree::SolveTT(int n) {
    double deltaT = maturity / n;
    double up = exp(sigma * sqrt(3 * deltaT));
    double pu = (1. / 6. + (r - q - pow(sigma, 2) / 2) * sqrt(deltaT / 12 / pow(sigma, 2))) * exp(-r * deltaT);
    double pm = 2. / 3. * exp(-r * deltaT);
    double pd = exp(-r * deltaT) - pu - pm;
    std::vector<double> price(2 * n + 1, 0);

    for (int i = 0; i <= 2 * n; i++) {
        if (!isCall)
            price[i] = K - S * pow(up, n - i);
        else
            price[i] = S * pow(up, n - i) - K;
        if (price[i] < 0)
            price[i] = 0;
    }

    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j <= 2 * i; j++) {
            price[j] = pu * price[j] + pm * price[j + 1] + pd * price[j + 2];
            if (!isCall && !isEuro) {
                double exercise = K - S * pow(up, i - j);
                if (price[j] < exercise)
                    price[j] = exercise;
            }
        }
        if (i == 1){
//            std::cout << price[0] << std::endl;
            delta = (price[0] - price[2]) / S / (up - 1 / up);
        }

        if (i == 2)
            gamma = ((price[0] - price[2]) / S / (up * up - 1)-(price[2] - price[4]) / S / (1 - 1 / up / up)) / S / (up - 1 / up);
    }
    theta = (price[1] - price[0]) / deltaT;
    return price[0];
}

double TrinoTree::SolveATT(int n) {
    return (SolveTT(n) + SolveTT(n - 1)) / 2;
}

double TrinoTree::SolveBTT(int n) {
    double deltaT = maturity / n;
    double up = exp(sigma * sqrt(3. * deltaT));
    double pu = (1. / 6. + (r - q - pow(sigma, 2) / 2.) * sqrt(deltaT / 12. / pow(sigma, 2.))) * exp(-r * deltaT);
    double pm = 2. / 3. * exp(-r * deltaT);
    double pd = exp(-r * deltaT) - pu - pm;
    std::vector<double> price(2 * n + 1, 0);

    for (int i = 0; i <= (2 * n - 2); i++) {
        BlackScholes temp(S * pow(up, n - 1 - i), K, deltaT, sigma, r, q);
        if (!isCall)
            price[i] = std::max(temp.PutPrice(), K - S * pow(up, n - 1 - i));
        else
            price[i] = temp.CallPrice();
    }

    for (int i = n - 2; i >= 0; i--) {
        for (int j = 0; j <= 2 * i; j++) {
            price[j] = (pu * price[j] + pm * price[j + 1] + pd * price[j + 2]);
            if (!isCall && !isEuro) {
                double exercise = K - S * pow(up, i - j);
                if (price[j] < exercise)
                    price[j] = exercise;
            }
        }
        if (i == 1){
//            std::cout << price[2] << std::endl;
            delta = (price[0] - price[2]) / S / (up - 1. / up);
        }

        if (i == 2)
            gamma = ((price[0] - price[2]) / S / (up * up - 1.)-(price[2] - price[4]) / S / (1. - 1. / up / up)) / S / (up - 1. / up);
    }
    theta = (price[1] - price[0]) / deltaT;
    return price[0];
}

double TrinoTree::SolveBTTR(int n) {
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

double TrinoTree::Delta(int n, int i) {
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

double TrinoTree::Gamma(int n, int i) {
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

double TrinoTree::Theta(int n, int i) {
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

TrinoTree::TrinoTree(const TrinoTree& orig) {
}

TrinoTree::~TrinoTree() {
}

