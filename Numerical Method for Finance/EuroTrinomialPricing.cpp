/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "EuroTrinomialPricing.h"
#include "math.h"
#include "BlackScholes.h"
#include <iostream>
#include <vector>
#include <cmath>
double V00,V12,V11, V10, V20, V22, V24, S00, S12,S11, S10, S24, S22, S20, delta, theta, Gamma;
EuroTrinoTree::EuroTrinoTree(double m, double s, double k,
        double sig, double Q, double R,
        bool call, bool euro) : maturity(m), S(s), K(k), sigma(sig), q(Q), r(R), isCall(call), isEuro(euro) {

}

void EuroTrinoTree::pickType(bool e, bool c) {
    isEuro = e;
    isCall = c;
}

std::vector<double> EuroTrinoTree::SolveTS(int n) {
	double deltaT = maturity / n;
	double up = exp(sigma * sqrt(3*deltaT));
	double down = 1.0 / up;
	double p0 = 1.0 / 6 + (r - q - sigma*sigma / 2.0)*sqrt(deltaT / (12 * sigma*sigma));
	double p1 = 2 / 3.0;
	double p2 = 1 / 3.0 - p0;
	std::vector<double> price(2*n + 1);

	for (int i = 0; i <= 2*n; i++) {
		if (!isCall)
			price[i] = K - S * pow(up, n-i);
		else
			price[i] = S * pow(up, n-i) - K;
		if (price[i] < 0)
			price[i] = 0;
	}

	for (int i = n - 1; i >= 0; i--) {
		for (int j = 0; j <= 2*i; j++) {
			price[j] = exp(-r*deltaT)*(p0 * price[j] + p1 * price[j + 1]+p2*price[j+2]);
			if (!isCall && !isEuro) {
				double exercise = K - S * pow(up, i - j);
				//                std::cout << price[j] - exercise << std::endl;
				if (price[j] < exercise)
					price[j] = exercise;
			}
		}
		if (i == 2) {
			V24 = price[4];
			V22 = price[2];
			V20 = price[0];
		}
		else if (i == 1) {
			V12 = price[2];
			V11 = price[1];
			V10 = price[0];
		}
		else if (i == 0) {
			V00 = price[0];
		}
	}
	S10 = up*S;
	S12 = down*S;
	S22 = S;
	S20 = up*up*S;
	S24 = down*down*S;
	delta = (V10 - V12) / (S10 - S12);
	Gamma = ((V20 - V22) / (S20 - S22) - (V22 - V24) / (S22 - S24)) / (S10 - S12);
	theta = (V11 - V00) / deltaT;
	return{ price[0],delta,Gamma,theta };
}


std::vector<double> EuroTrinoTree::SolveTBS(int n) {
	double deltaT = maturity / n;
	double up = exp(sigma * sqrt(3 * deltaT));
	double down = 1.0 / up;
	double p0 = 1.0 / 6 + (r - q - sigma*sigma / 2.0)*sqrt(deltaT / (12 * sigma*sigma));
	double p1 = 2 / 3.0;
	double p2 = 1 / 3.0 - p0;
	std::vector<double> price(2*n + 1);


	for (int i = 0; i <= 2*n - 2; i++) {
		BlackScholes temp(S * pow(up,n-i-1), K, deltaT, sigma, r, q);
		if (!isCall)
			price[i] = temp.PutPrice();
		else
			price[i] = temp.CallPrice();
	}

	for (int i = n - 2; i >= 0; i--) {
		for (int j = 0; j <= 2*i; j++) {
			price[j] = exp(-r*deltaT)*(p0 * price[j] + p1 * price[j + 1] + p2*price[j + 2]);
			if (!isCall && !isEuro) {
				double exercise = K - S * pow(up, i - j);
				//                std::cout << price[j] - exercise << std::endl;
				if (price[j] < exercise)
					price[j] = exercise;
			}
		}
		if (i == 2) {
			V24 = price[4];
			V22 = price[2];
			V20 = price[0];
		}
		else if (i == 1) {
			V12 = price[2];
			V11 = price[1];
			V10 = price[0];
		}
		else if (i == 0) {
			V00 = price[0];
		}
	}
	S10 = up*S;
	S12 = down*S;
	S22 = S;
	S20 = up*up*S;
	S24 = down*down*S;
	delta = (V10 - V12) / (S10 - S12);
	Gamma = ((V20 - V22) / (S20 - S22) - (V22 - V24) / (S22 - S24)) / (S10 - S12);
	theta = (V11 - V00) / deltaT;
	return{ price[0],delta,Gamma,theta };
}

std::vector<double> EuroTrinoTree::SolveTBSR(int n) {
	auto vec1 = SolveTBS(n);
	auto vec2 = SolveTBS(n/2);
	std::vector<double> result;
	for (int i = 0; i < 4; ++i)
		result.push_back(2*vec1[i]-vec2[i]);
	return result;
}

