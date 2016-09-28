/*
 * Functions.h
 *
 *  Created on: Aug 27, 2016
 *      Author: wenhaohu
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <Eigen/Dense>
#include <utility>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <tuple>

using Eigen::internal::BandMatrix;
using Eigen::Dynamic;

Eigen::VectorXd forward_subst(Eigen::MatrixXd, Eigen::VectorXd);
Eigen::VectorXd forward_subst(Eigen::MatrixXd, Eigen::VectorXd, int band);
Eigen::VectorXd backward_subst(Eigen::MatrixXd, Eigen::VectorXd);
Eigen::VectorXd backward_subst(Eigen::MatrixXd, Eigen::VectorXd, int band);
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> lu_no_pivoting(Eigen::MatrixXd);
std::tuple<Eigen::PermutationMatrix<Dynamic, Dynamic>, Eigen::MatrixXd, Eigen::MatrixXd> lu_row_pivoting(Eigen::MatrixXd);
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> lu_no_pivoting_band(Eigen::MatrixXd,
        int);
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> lu_row_pivoting_band(
        Eigen::MatrixXd, int);
Eigen::MatrixXd cholesky_spd(Eigen::MatrixXd);
Eigen::MatrixXd cholesky_band_spd(Eigen::MatrixXd, int);
Eigen::MatrixXd ReadMatrix(std::string);


#endif /* FUNCTIONS_H_ */
