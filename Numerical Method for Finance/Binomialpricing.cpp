/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "BinomialPricing.h"
#include "math.h"
#include "BlackScholes.h"
#include <iostream>
#include <vector>

BinoTree::BinoTree(double m, double s, double k,
        double sig, double Q, double R,
        bool call, bool euro) : maturity(m), S(s), K(k), sigma(sig), q(Q), r(R), isCall(call), isEuro(euro) {

}

void BinoTree::pickType(bool e, bool c) {
    isEuro = e;
    isCall = c;
}

double BinoTree::SolveBS(int n) {
    double deltaT = maturity / n;
    double up = exp(sigma * sqrt(deltaT));
    double p0 = (up * exp(-r * deltaT) - exp(-q * deltaT)) * up / (pow(up, 2) - 1);
    double p1 = exp(-r * deltaT) - p0;
    std::vector<double> price(n + 1, 0);

    for (int i = 0; i <= n; i++) {
        if (!isCall)
            price[i] = K - S * pow(up, 2 * i - n);
        else
            price[i] = S * pow(up, 2 * i - n) - K;
        if (price[i] < 0)
            price[i] = 0;
    }

    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j <= i; j++) {
            price[j] = p0 * price[j] + p1 * price[j + 1];
            if (!isCall && !isEuro) {
                double exercise = K - S * pow(up, 2 * j - i);
                if (price[j] < exercise)
                    price[j] = exercise;
            }
        }
    }
    return price[0];
}

double BinoTree::SolveABS(int n) {
    return (SolveBS(n) + SolveBS(n - 1)) / 2;
}

double BinoTree::SolveBBS(int n) {
    double deltaT = maturity / n;
    double up = exp(sigma * sqrt(deltaT));
    double p0 = (up * exp(-r * deltaT) - exp(-q * deltaT)) * up / (pow(up, 2) - 1);
    double p1 = exp(-r * deltaT) - p0;
    std::vector<double> price(n + 1, 0);

    for (int i = 0; i <= n - 1; i++) {
        BlackScholes temp(S * pow(up, 2 * i - n + 1), K, deltaT, sigma, r, q);
        if (!isCall)
            price[i] = temp.PutPrice();
        else
            price[i] = temp.CallPrice();
    }

    for (int i = n - 2; i >= 0; i--) {
        for (int j = 0; j <= i; j++) {
            price[j] = p0 * price[j] + p1 * price[j + 1];
            if (!isCall && !isEuro) {
                double exercise = K - S * pow(up, 2 * j - i);
                if (price[j] < exercise)
                    price[j] = exercise;
            }
        }
    }
    return price[0];
}

double BinoTree::SolveBBSR(int n) {
    return 2 * SolveBBS(n) - SolveBBS(n / 2);
}