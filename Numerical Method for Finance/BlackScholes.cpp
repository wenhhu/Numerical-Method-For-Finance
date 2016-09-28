/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BlackScholes.cpp
 * Author: wenhaohu
 * 
 * Created on January 3, 2016, 10:03 PM
 */

#include "BlackScholes.h"
#include "NumMethod.h"
#include "math.h"

BlackScholes::BlackScholes(double s, double k, double t, double ssigma, double R, double Q)\
:S(s),K(k),T(t),sigma(ssigma),r(R),q(Q) {
}

BlackScholes::BlackScholes(const BlackScholes& orig) {
}

BlackScholes::~BlackScholes() {
}

double BlackScholes::CallPrice(double t){
    double d1=(log(S/K)+(r-q+sigma*sigma/2)*(T-t))/sigma/sqrt(T-t);
    double d2=d1-sigma*sqrt(T-t);
    return S*exp(-q*(T-t))*cum_dist_normal(d1)-K*exp(-r*(T-t))*cum_dist_normal(d2);
}

double BlackScholes::PutPrice(double t){
    double d1=(log(S/K)+(r-q+sigma*sigma/2)*(T-t))/sigma/sqrt(T-t);
    double d2=d1-sigma*sqrt(T-t);
    return -S*exp(-q*(T-t))*cum_dist_normal(-d1)+K*exp(-r*(T-t))*cum_dist_normal(-d2);
}

BlackScholes::BlackScholes(){
    
}