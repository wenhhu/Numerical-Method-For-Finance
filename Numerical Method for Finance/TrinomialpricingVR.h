/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TrinomialpricingVR.h
 * Author: wenhaohu
 *
 * Created on October 12, 2016, 11:29 PM
 */

#ifndef TRINOMIALPRICINGVR_H
#define TRINOMIALPRICINGVR_H

class TrinoTreeVR{
    double maturity;
    double S, K, sigma, q, r, delta, gamma, theta;
    bool isCall, isEuro;
    
public:
    //double m, double s, double k, double sig, double Q, double R, bool call, bool euro
    TrinoTreeVR(double, double, double, double, double, double, bool, bool);
    void pickType(bool, bool);
    double SolveTT(int);
    double SolveATT(int);
    double SolveBTT(int);
    double SolveBTTR(int);
    double Delta(int,int);
    double Gamma(int,int);
    double Theta(int,int);
    TrinoTreeVR();
    TrinoTreeVR(const TrinoTreeVR& orig);
    ~TrinoTreeVR();


};


#endif /* TRINOMIALPRICINGVR_H */

