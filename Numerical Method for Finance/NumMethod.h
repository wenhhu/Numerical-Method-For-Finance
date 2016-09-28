/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   NumMethod.h
 * Author: wenhaohu
 *
 * Created on January 1, 2016, 6:51 PM
 */

#ifndef NUMMETHOD_H
#define NUMMETHOD_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif



#define PI 3.141592653

//Integration algorithm
double TrapezoidalRule(double (*f)(double), double error, double a, double b);
double MidpointRule(double (*f)(double), double error, double a, double b);
double SimpsonRule(double (*f)(double), double error, double a, double b);

//Bond
double BondP_ZRC(int n, double* t, double* cash, double (*f)(double));
double BondP_IIRC(int n, double* t, double* cash, double (*f)(double), double* tol);
double BondP_Yield(int n, double* t, double* cash, double y);

//Cumulative distribution
double cum_dist_normal(double t);
double cum_dist_normal_simp(double t);

//Root finding
double NewtonMethod(double (*f)(double), double x0, double tol_approx, double tol_consec);
double BisectionMethod(double (*f)(double), double left, double right, double tol_approx, double tol_int);

#endif /* NUMMETHOD_H */

