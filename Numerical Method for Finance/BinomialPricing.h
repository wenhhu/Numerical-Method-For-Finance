/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BinomialPricing.h
 * Author: wenhaohu
 *
 * Created on September 17, 2016, 12:03 PM
 */

#ifndef BINOMIALPRICING_H
#define BINOMIALPRICING_H

class BinoTree{
    double maturity;
    double S, K, sigma, q, r;
    bool isCall, isEuro;


public:
//sdouble m, double s, double k, double sig, double Q, double R, bool call, bool euro
    BinoTree(double, double, double, double, double, double, bool, bool);
    void pickType(bool, bool);
    double SolveBS(int);
    double SolveABS(int);
    double SolveBBS(int);
    double SolveBBSR(int);
};



#endif /* BINOMIALPRICING_H */

