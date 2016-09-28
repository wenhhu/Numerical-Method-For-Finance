/*
 * Functions.cpp
 *
 *  Created on: Aug 27, 2016
 *      Author: wenhaohu
 */

#include "Functions.h"
#include <algorithm>

using Eigen::internal::BandMatrix;
using std::max;
using std::min;
using Eigen::Dynamic;
using Eigen::PermutationMatrix;

Eigen::VectorXd forward_subst(Eigen::MatrixXd L, Eigen::VectorXd b) {
    assert(L.rows() == L.cols() && L.rows() == b.rows());
    int n = L.rows();
    Eigen::VectorXd res(n);
    double temp;
    res(0) = b(0) / L(0, 0);
    for (int i = 1; i < n; i++) {
        temp = 0;
        for (int j = 0; j < i; j++)
            temp += L(i, j) * res(j);
        res(i) = (b(i) - temp) / L(i, i);
    }
    return res;
}

Eigen::VectorXd forward_subst(Eigen::MatrixXd L, Eigen::VectorXd b, int m) {
    assert(L.rows() == L.cols() && L.rows() == b.rows());
    int n = L.rows();
    Eigen::VectorXd res(n);
    double temp;
    res(0) = b(0) / L(0, 0);
    for (int i = 1; i < n; i++) {
        temp = 0;
        for (int j = max(i - m, 0); j < i; j++)
            temp += L(i, j) * res(j);
        res(i) = (b(i) - temp) / L(i, i);
    }
    return res;
}

Eigen::VectorXd backward_subst(Eigen::MatrixXd U, Eigen::VectorXd b) {
    assert(U.rows() == U.cols() && U.rows() == b.rows());
    int n = U.rows();
    Eigen::VectorXd res(n);
    double temp;
    res(n - 1) = b(n - 1) / U(n - 1, n - 1);
    for (int i = n - 2; i >= 0; i--) {
        temp = 0;
        for (int j = n - 1; j > i; j--)
            temp += U(i, j) * res(j);
        res(i) = (b(i) - temp) / U(i, i);
    }
    return res;
}

Eigen::VectorXd backward_subst(Eigen::MatrixXd U, Eigen::VectorXd b, int m) {
    assert(U.rows() == U.cols() && U.rows() == b.rows());
    int n = U.rows();
    Eigen::VectorXd res(n);
    double temp;
    res(n - 1) = b(n - 1) / U(n - 1, n - 1);
    for (int i = n - 2; i >= 0; i--) {
        temp = 0;
        for (int j = min(n - 1, i + m); j > i; j--)
            temp += U(i, j) * res(j);
        res(i) = (b(i) - temp) / U(i, i);
    }
    return res;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> lu_no_pivoting(Eigen::MatrixXd A) {
    assert(A.rows() == A.cols());
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> res;
    Eigen::MatrixXd L(A.rows(), A.cols()), U(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = i; j < A.rows(); j++) {
            U(i, j) = A(i, j);
            L(j, i) = A(j, i) / U(i, i);
        }
        for (int k = i + 1; k < A.rows(); k++) {
            for (int l = i + 1; l < A.rows(); l++)
                A(k, l) = A(k, l) - L(k, i) * U(i, l);
        }
    }
    res.first = L;
    res.second = U;
    return res;
}

std::tuple<PermutationMatrix<Dynamic, Dynamic>, Eigen::MatrixXd, Eigen::MatrixXd> lu_row_pivoting(Eigen::MatrixXd A) {
    assert(A.rows() == A.cols());
    std::tuple<PermutationMatrix<Dynamic, Dynamic>, Eigen::MatrixXd, Eigen::MatrixXd> res;
    Eigen::MatrixXd L(A.rows(), A.cols()), U(A.rows(), A.cols());
    PermutationMatrix<Dynamic, Dynamic> P(A.rows());
    P.setIdentity();
    int max_col = 0;
    for (int i = 0; i < A.rows(); i++) {
        //check the index of row with maximum A(j,i)
        max_col = i;
        for (int j = i; j < A.rows(); j++) {
            if (fabs(A(max_col, i)) < fabs(A(j, i)))
                max_col = j;
        }
        P.applyTranspositionOnTheLeft(i, max_col);
        PermutationMatrix<Dynamic, Dynamic> temp(A.rows());
        temp.setIdentity();
        temp.applyTranspositionOnTheLeft(i, max_col);
        //pivoting of L
        L = temp * L;

        //pivoting of A
        A = temp * A;

        //LU decomposition
        for (int j = i; j < A.rows(); j++) {
            U(i, j) = A(i, j);
            L(j, i) = A(j, i) / U(i, i);
        }
        for (int k = i + 1; k < A.rows(); k++) {
            for (int l = i + 1; l < A.rows(); l++)
                A(k, l) = A(k, l) - L(k, i) * U(i, l);
        }
    }
    return std::make_tuple(P, L, U);
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> lu_no_pivoting_band(Eigen::MatrixXd A, int m) {
    assert(A.rows() == A.cols());
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> res;
    Eigen::MatrixXd L(A.rows(), A.cols()), U(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = i; j < min(int(A.rows()), i + m + 1);
                j++) {
            U(i, j) = A(i, j);
            L(j, i) = A(j, i) / U(i, i);
        }
        for (int k = i + 1; k < A.rows(); k++) {
            for (int l = max(i + 1, k - m);
                    l < min(int(A.rows()), k + m + 1);
                    l++) {
                A(k, l) = A(k, l) - L(k, i) * U(i, l);
            }
        }
    }
    res.first = L;
    res.second = U;
    return res;
}

Eigen::MatrixXd ReadMatrix(std::string name) {
    std::string line, val;
    std::stringstream iss;
    std::ifstream file(name.data());
    int m = 0, n = 0;
    while (file.good()) {
        getline(file, line);
        iss << line;
        while (getline(iss, val, ','))
            n++;
        m++;
    }
    file = std::ifstream(name.data());
    Eigen::MatrixXd res(m, n);
    for (int i = 0; file.good(); i++) {
        getline(file, line);
        iss = std::stringstream(line);
        for (int j = 0; iss.good(); j++) {
            getline(iss, val, ',');
            if (!val.empty() & int(val[0]) != 13) {
                res(i, j) = std::stod(val);
            } else
                res(i, j) = 0;
        }
    }
    return res;
}

Eigen::MatrixXd cholesky_spd(Eigen::MatrixXd A) {
    assert(A.transpose() == A);
    Eigen::MatrixXd U(A.rows(), A.cols());
    for (int i = 0; i < A.rows() - 1; i++) {
        U(i, i) = sqrt(A(i, i));
        for (int j = i; j < A.rows(); j++)
            U(i, j) = A(i, j) / U(i, i);
        for (int k = i + 1; k < A.rows(); k++) {
            for (int l = i + 1; l < A.rows(); l++)
                A(k, l) = A(k, l) - U(i, k) * U(i, l);
        }
    }
    U(A.rows() - 1, A.rows() - 1) = sqrt(A(A.rows() - 1, A.rows() - 1));
    return U;
}

Eigen::MatrixXd cholesky_band_spd(Eigen::MatrixXd A, int m) {
    assert(A.transpose() == A);
    Eigen::MatrixXd U(A.rows(), A.cols());
    for (int i = 0; i < A.rows() - 1; i++) {
        U(i, i) = sqrt(A(i, i));
        for (int j = i; j < min(i + m + 1, int(A.rows())); j++)
            U(i, j) = A(i, j) / U(i, i);
        for (int k = i + 1; k < A.rows(); k++) {
            for (int l = max(i + 1, k - m); l < min(int(A.rows()), k + m + 1); l++)
                A(k, l) = A(k, l) - U(i, k) * U(i, l);
        }
    }
    U(A.rows() - 1, A.rows() - 1) = sqrt(A(A.rows() - 1, A.rows() - 1));
    return U;
}
