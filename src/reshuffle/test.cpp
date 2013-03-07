/*
 * test.cpp
 *
 *  Created on: 07.03.2013
 *      Author: lima
 */
#include "test.h"
#include <string.h>
#include <stdlib.h>
using namespace std;

parameters_test::parameters_test(int a){
	cout<<"START TEST"<<endl;
	string cmdline_1="4test--datadims--snps=1-100,50-150,229--traits=1-3,2-5,10";
	int argc_1=2;
	char *argv_1[2];
	argv_1[1] = new char[cmdline_1.length()+1];
	strcpy(argv_1[1], cmdline_1.c_str());
	Params_1=Parameters(argc_1,argv_1);
}


