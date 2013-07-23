/*
 * test.h
 *
 *  Created on: 07.03.2013
 *      Author: lima
 */

#ifndef TEST_H_
#define TEST_H_

typedef long long int64_t;


#include "Parameters.h"
#include "iout_file.h"
#include "reshuffle.h"
//#include <dir.h>
using namespace std;

class test{
public:
	string name;
	string result;
	string check;
	Parameters *Params_test;//&Paramss;
	test(string,string,string,string);
	void run(ofstream&);
};
#endif
/* TEST_H_ */
