/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TrinomialPricing.h
 * Author: Zhihao Chen
 *
 * Created on October 13, 2016, 12:03 PM
 */

#ifndef EUROTRINOMIALPRICING_H
#define EUROTRINOMIALPRICING_H
#include <vector>
class EuroTrinoTree{
    double maturity;
    double S, K, sigma, q, r;
    bool isCall, isEuro;


public:
//sdouble m, double s, double k, double sig, double Q, double R, bool call, bool euro
    EuroTrinoTree(double, double, double, double, double, double, bool, bool);
    void pickType(bool, bool);
    std::vector<double> SolveTS(int);
    std::vector<double> SolveTBS(int);
    std::vector<double> SolveTBSR(int);
};



#endif /* EUROTRINOMIALPRICING_H */

