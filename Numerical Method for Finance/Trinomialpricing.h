/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Trinomialpricing.h
 * Author: wenhaohu
 *
 * Created on October 12, 2016, 4:48 PM
 */

#ifndef TRINOMIALPRICING_H
#define TRINOMIALPRICING_H

#include "BinomialPricing.h"

class TrinoTree{
    double maturity;
    double S, K, sigma, q, r, delta, gamma, theta;
    bool isCall, isEuro;
    
public:
    //double m, double s, double k, double sig, double Q, double R, bool call, bool euro
    TrinoTree(double, double, double, double, double, double, bool, bool);
    void pickType(bool, bool);
    double SolveTT(int);
    double SolveATT(int);
    double SolveBTT(int);
    double SolveBTTR(int);
    double Delta(int,int);
    double Gamma(int,int);
    double Theta(int,int);
    TrinoTree();
    TrinoTree(const TrinoTree& orig);
    ~TrinoTree();


};

#endif /* TRINOMIALPRICING_H */

