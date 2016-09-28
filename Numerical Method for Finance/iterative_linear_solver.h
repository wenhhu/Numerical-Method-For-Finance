/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   iterative_linear_solver.h
 * Author: wenhaohu
 *
 * Created on September 11, 2016, 9:29 AM
 */

#ifndef ITERATIVE_LINEAR_SOLVER_H
#define ITERATIVE_LINEAR_SOLVER_H

#include <Eigen/Dense>

class IterativeSolver {
protected:
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    bool isResStop;
    Eigen::VectorXd x0;
    int ic = 0;
    double tol;
public:
    IterativeSolver(Eigen::MatrixXd, Eigen::VectorXd, double, bool pick = true);
    IterativeSolver() = default;
    void setA(Eigen::MatrixXd);
    void setb(Eigen::VectorXd);
    int getic();
    void setStop(bool);
    virtual Eigen::VectorXd Solve(Eigen::VectorXd) = 0;
};

class Jacobi:public IterativeSolver{
    using IterativeSolver::IterativeSolver;
public:
    virtual Eigen::VectorXd Solve(Eigen::VectorXd);
};

class GaussSeidel:public IterativeSolver{
    using IterativeSolver::IterativeSolver;
public:
    Eigen::VectorXd getA();
    virtual Eigen::VectorXd Solve(Eigen::VectorXd);
};

class SOR:public IterativeSolver{
    double omega;
public:
    SOR(Eigen::MatrixXd, Eigen::VectorXd, double, double, bool pick = true);
    void setOmega(double);
    double getOmega();
    Eigen::MatrixXd getR();
    virtual Eigen::VectorXd Solve(Eigen::VectorXd);
};


#endif /* ITERATIVE_LINEAR_SOLVER_H */

