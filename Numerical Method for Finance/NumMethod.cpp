/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "NumMethod.h"
#include <iostream>
#include <Eigen/Dense>


double TrapezoidalRule(double (*f)(double), double error, double a, double b){
    int n=4;
    double h = (b-a)/n;
    double res = ((*f)(a)+(*f)(b))/2;
    double temp = res;
    
//    #pragma omp parallel for reduction(+: res, temp)
    for(int i=1; i<n; i++){
        res += (*f)(a+i*h);
        temp += (*f)(a+(2*i-1)*h/2);
        temp += (*f)(a+2*i*h/2);
    }
    
    temp += (*f)(a+(2*n-1)*h/2);
    res *= h;
    temp = temp*h/2;

    double dif = fabs(temp-res);
    while(dif>=error){
        n *= 2;
        h /=2;
        res = temp;
        temp = ((*f)(a)+(*f)(b))/2;
    
        for(int i=1; i<n; i++){
            temp += (*f)(a+(2*i-1)*h/2);
            temp += (*f)(a+2*i*h/2);
        }
        
        temp += (*f)(a+(2*n-1)*h/2);
        temp = temp*h/2;
        
//        std::cout << "n=" << n/2 << ": " << res << std::endl;
        dif = fabs(temp-res);
    }

    return temp;
}

double MidpointRule(double (*f)(double), double error, double a, double b){
    int n = 1;
    double h = (b-a)/n;
    double res = 0;
    double temp = 0;
    bool flag=1;
    for(int i=0;i<n;i++)
        res+=(f(a+h/2+i*h)*h);
    while(flag){
        n *= 2;
        h = (b-a)/n;
        temp = res;
        res = 0;
        for(int i=0;i<n;i++)
            res+=(f(a+h/2+i*h)*h);
        flag = (fabs(res - temp) >= error);
    }
    return res;
}

double SimpsonRule(double (*f)(double), double error, double a, double b){
    int n = 4;
    double h = (b-a)/n;
    double res = 0, even = 0, odd = 0;
    double temp = 0;
    int iter=1;
    bool flag=1;
    for(int i=0; i < n/2; i++){
        odd += 4*f(a+(2*i+1)*h);
        even += 2*f(a+2*(i+1)*h);
    }
    res = (f(a)+odd+even+f(b))*h/3;
//    std::cout.precision(12);
//    std::cout << iter << ": " << '\n' << std::fixed << 0.5+1/sqrt(2*PI)*res << std::endl;
    while(flag){
        iter++;
        n *= 2;
        h = (b-a)/n;
        temp = res;
        res = 0;
        for(int i=0;i<n;i++)
            res+=(f(a+h/2+i*h)*h);
//        std::cout << iter << ": " << '\n' << std::fixed << 0.5+1/sqrt(2*PI)*res << std::endl;
        flag = (fabs(res - temp) >= error);
    }
    return res;
}

double BondP_ZRC(int n, double* t, double* cash, double (*r)(double)){
    double res=0;
    double dist=0;
    for(int i=0;i<n;i++){
        dist = exp(-r(t[i])*t[i]);
        res += dist*cash[i];
    }
    return res;
}

double BondP_IIRC(int n, double* t, double* cash, double (*f)(double), double* tol){
    double res=0;
    double dist=0;
    for(int i=0;i<n;i++){
        dist = exp(-SimpsonRule(f,tol[i],0,t[i]));
        res += dist*cash[i];
    }
    return res;
}

double BondP_Yield(int n, double* t, double* cash, double y){
    double res[3]={0,0,0};
    double dist=0;
    for(int i=0;i<n;i++){
        dist = exp(-t[i]*y);
        res[0] += dist*cash[i];
        res[1] += dist*cash[i]*t[i];
        res[2] += dist*cash[i]*t[i]*t[i];
    }
    res[1] /= res[0];
    res[2] /= res[0];
    return res[0];
}

double cum_dist_normal(double t){
    double z=fabs(t);
    double y=1/(1+0.2316419*z);
    double a1 = 0.319381530, a2 = -0.356563782, a3 = 1.781477937, a4 = -1.821255978, a5 = 1.330274429;
    double m=1-exp(-t*t/2)*(a1*y+a2*y*y+a3*pow(y,3)+a4*pow(y,4)+a5*pow(y,5))/sqrt(2*PI);
    if(t>0)
        return m;
    else
        return 1-m;
}

double pdf(double x){
    return exp(-x*x/2);
}

double cum_dist_normal_simp(double t){
    double res=0.5;
    if(t>0)
        res += 1/sqrt(2*PI)*SimpsonRule(pdf,1e-12,0,t);
    else
        res += 1/sqrt(2*PI)*SimpsonRule(pdf,1e-12,t,0);
    return res;
}

double NewtonMethod(double (*f)(double), double x0, double tol_approx, double tol_consec){
    double xnew=x0;
    double xold=x0-1;
    std::cout << xnew << std::endl;
    while(fabs(f(xnew))>tol_approx || fabs(xnew-xold)>tol_consec){
        xold=xnew;
        xnew=xold-f(xold)/(f(xold+tol_consec)-f(xold-tol_consec))*2*tol_consec;
        std::cout << xnew << std::endl;
    }
    return xnew;
}

double BisectionMethod(double (*f)(double), double left, double right, double tol_approx, double tol_int){
    double xl=left,xr=right;
    double xm;
    while(std::max(fabs(f(xl)),fabs(f(xr)))>tol_approx || fabs(xr-xl)>tol_int){
        xm=(xl+xr)/2;
        if(f(xl)*f(xm)<0)
            xr=xm;
        else
            xl=xm;
    }
    return xl;
}