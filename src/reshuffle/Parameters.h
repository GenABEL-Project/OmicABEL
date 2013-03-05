/*
 * Parameters.h
 *
 *  Created on: 28.02.2013
 *      Author: lima
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <string>
#include <set>
using namespace std;

class Parameter {

public:
	string name; 	//name of parametr
	bool use; 	// Parameter use or not (Is paameter in command line?)
	string value; 	// value of parametr,chars after "=" symbol
	Parameter(string,string); 	//constructor
	Parameter();		//default constructor
	set<int> valueset;
};

ostream &operator <<(ostream &, Parameter);

class Parameters {
public:
	string iout_fname;
	string out_fname;
	Parameter datadims;
	Parameter snpnames;
	Parameter traitnames;
	Parameter traits;
	Parameter snps;
	Parameter heritabilities;
	Parameter chi;
	Parameter dataslim;
	Parameters(int, char*[]);		//	Constructor from cmdline
	string get_cmd_line(int,char*[]);

};

ostream &operator <<(ostream &, Parameters);

#endif /* PARAMETERS_H_ */
