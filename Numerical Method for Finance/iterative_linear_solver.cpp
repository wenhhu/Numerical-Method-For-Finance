/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "iterative_linear_solver.h"
#include "Functions.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

IterativeSolver::IterativeSolver(MatrixXd Ai,
        VectorXd bi,
        double t,
        bool pick) : A(Ai), b(bi), isResStop(pick), tol(t) {

}

void IterativeSolver::setA(MatrixXd Ai) {
    A = Ai;
}

void IterativeSolver::setb(VectorXd bi) {
    b = bi;
}

void IterativeSolver::setStop(bool pick) {
    isResStop = pick;
}

int IterativeSolver::getic() {
    return ic;
}

VectorXd Jacobi::Solve(VectorXd init) {
    ic = 0;
    x0 = init;
    auto x = x0;
    VectorXd r = b - A*x0;
    double stop_iter_resid = tol * r.norm();
    auto D = MatrixXd((A.diagonal()).asDiagonal());
    auto L = MatrixXd(A.triangularView<Eigen::StrictlyLower>());
    auto U = MatrixXd(A.triangularView<Eigen::StrictlyUpper>());
    VectorXd b_new = forward_subst(D, b, 0);

    if (isResStop) {
        while (r.norm() > stop_iter_resid) {
            x = -forward_subst(D, (L + U) * x, 0) + b_new;
            r = b - A*x;
            if(ic < 3){
                std::cout << "The " << ic+1 << "th iteration is:" << std::endl;
                std::cout << x << std::endl << std::endl;
            }
            ic++;
        }
    } else {
        VectorXd temp = VectorXd::Zero((A.rows()));
        while ((x - temp).norm() > tol) {
            temp = x;
            x = -forward_subst(D, (L + U) * x, 0) + b_new;
            ic++;
        }
    }
    return x;
}

VectorXd GaussSeidel::Solve(VectorXd init) {
    ic = 0;
    x0 = init;
    auto x = x0;
    VectorXd r = b - A*x0;
    auto stop_iter_resid = tol * r.norm();
    auto D = MatrixXd((A.diagonal()).asDiagonal());
    auto L = MatrixXd(A.triangularView<Eigen::StrictlyLower>());
    auto U = MatrixXd(A.triangularView<Eigen::StrictlyUpper>());
    VectorXd b_new = forward_subst(D + L, b);

    if (isResStop) {
        while (r.norm() > stop_iter_resid) {
            x = -forward_subst(D + L, U * x) + b_new;
            r = b - A*x;
            if(ic < 3){
                std::cout << "The " << ic+1 << "th iteration is:" << std::endl;
                std::cout << x << std::endl << std::endl;
            }
            ic++;
        }
    } else {
        VectorXd temp = VectorXd::Zero((A.rows()));
        while ((x - temp).norm() > tol) {
            temp = x;
            x = -forward_subst(D + L, U * x) + b_new;
            ic++;
        }
    }
    return x;
}

Eigen::VectorXd GaussSeidel::getA(){
    return A;
}

SOR::SOR(MatrixXd Ai,
        VectorXd bi,
        double t,
        double k,
        bool pick) : IterativeSolver(Ai, bi, t, pick), omega(k) {
}

Eigen::VectorXd SOR::Solve(Eigen::VectorXd init) {
    ic = 0;
    x0 = init;
    auto x = x0;
    VectorXd r = b - A*x0;
    auto stop_iter_resid = tol * r.norm();
    auto D = MatrixXd((A.diagonal()).asDiagonal());
    auto L = MatrixXd(A.triangularView<Eigen::StrictlyLower>());
    auto U = MatrixXd(A.triangularView<Eigen::StrictlyUpper>());
    VectorXd b_new = omega * forward_subst(D + omega*L, b);

    if (isResStop) {
        while (r.norm() > stop_iter_resid) {
            x = forward_subst(D + omega*L, (1 - omega) * D * x - omega * U * x) + b_new;
            r = b - A*x;
            if(ic < 3){
                std::cout << "The " << ic+1 << "th iteration is:" << std::endl;
                std::cout << x << std::endl << std::endl;
            }
            ic++;
        }
    } else {
        VectorXd temp = VectorXd::Zero((A.rows()));
        while ((x - temp).norm() > tol) {
            temp = x;
            x = forward_subst(D + omega*L, (1 - omega) * D * x - omega * U * x) + b_new;
            ic++;
        }
    }
    return x;
}

double SOR::getOmega() {
    return omega;
}

MatrixXd SOR::getR() {
    int n = A.rows();
    VectorXd x = VectorXd::Ones(3);
    auto D = MatrixXd((A.diagonal()).asDiagonal());
    auto L = MatrixXd(A.triangularView<Eigen::StrictlyLower>());
    auto U = MatrixXd(A.triangularView<Eigen::StrictlyUpper>());
    VectorXd b_new = omega * forward_subst(D + omega*L, b);
    return b_new;
}

void SOR::setOmega(double o){
    omega = o;
}