/*
 * main.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Содбо
 */
#include <iostream>
#include <ctime>
#include "iout_file.h"
#include "Parameters.h"
#include "Reshuffle.h"
#include <iterator>


using namespace std;

int main(int argc, char* argv[]) {
	cout << "Every day I'm [re]shuffling!" << endl;
	Parameters Params(argc, argv);
	cout << Params;
	iout_file iout_F(Params);
	cout << "finish iout_file read " << double(clock()) / CLOCKS_PER_SEC << endl;
	cout<<iout_F.header;
	cout<<iout_F.labels;
	Reshuffle reshh(iout_F,Params);
	reshh.run();
	cout << "finish_reshuffling " << double(clock()) / CLOCKS_PER_SEC << endl;
	return 0;
}
