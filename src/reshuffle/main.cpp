/*
 * main.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Sodbo
 */
#include <iostream>
#include <ctime>
#include "iout_file.h"
#include "parameters.h"
#include "Reshuffle.h"
#include <iterator>
#include "test.h"
//TODO Create Makefile

using namespace std;

int main(int argc, char* argv[]) {
	cout << "Every day I'm [re]shuf1fling!" << endl;
	parameters_test test(2);
	Parameters Params(argc, argv);
	if(Params.info.use)
		cout << Params;
	cout<<"SHETT";
	iout_file iout_F(Params);
	cout << "finish iout_file read " << double(clock()) / CLOCKS_PER_SEC << endl;
	if(Params.info.use){
		cout<<iout_F.header;
		cout<<iout_F.labels;
	}
	Reshuffle reshh(iout_F,Params);
	reshh.run();
	cout << "finish_reshuffling " << double(clock()) / CLOCKS_PER_SEC << endl;
	return 0;
}
