/*
 * Cholesky.h
 *
 *  Created on: Sep 4, 2016
 *      Author: wenhaohu
 */

#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#include "Functions.h"

class Cholesky {
	int m;
	Eigen::MatrixXd input;
	bool isBand;
	bool isTrid;

public:
	Cholesky();
	virtual ~Cholesky();
};

#endif /* CHOLESKY_H_ */
