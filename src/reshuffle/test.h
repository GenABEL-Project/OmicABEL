/*
 * test.h
 *
 *  Created on: 07.03.2013
 *      Author: lima
 */

#ifndef TEST_H_
#define TEST_H_
#include "parameters.h"

using namespace std;

class parameters_test{

public:
	Parameters Params_1;
	Parameters Params_2;
	Parameters Params_3;
	parameters_test(int);
	void run();
};

void run();
#endif /* TEST_H_ */
