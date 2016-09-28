/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BlackScholes.h
 * Author: wenhaohu
 *
 * Created on January 3, 2016, 10:03 PM
 */

#ifndef BLACKSCHOLES_H
#define BLACKSCHOLES_H


class BlackScholes {
public:
    BlackScholes(double S, double K, double T, double sigma, double r, double q);
    BlackScholes();
    BlackScholes(const BlackScholes& orig);
    double CallPrice(double t=0);
    double PutPrice(double t=0);
    virtual ~BlackScholes();

protected:
    double test;

private:
    double S, K, T, sigma, r, q;

};




#endif /* BLACKSCHOLES_H */

