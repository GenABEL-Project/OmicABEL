/*
 * Parameters.h
 *
 *  Created on: 28.02.2013
 *      Author: Sodbo
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <getopt.h>
#include <stdlib.h>

using namespace std;

class Parameter {

public:
	string name; 	//name of parameter
	bool use; 	// Parameter use or not (Is parameter in command line?)
	string value; 	// value of parameter,chars after "=" symbol
	set<int> numbersset;
	set<string> namesset;
	string delfromcmdline(string);
	void setbynames(vector<string>);
	string outfile;
	Parameter(string,string,string); 	//constructor
	Parameter();		//default constructor
};

ostream &operator <<(ostream &, Parameter);

class Parameters {
public:
	string iout_fname; // OmicABEL iout_file_name
	string out_fname; // OmicABEL out_file_name
	string outfile; // out_file path
	bool default_outfile; // use default out_file_name
	bool get_help; //Print help
	bool get_info; // Write info about programm's run
	bool run_test;
	bool write_datadims;
	bool write_slim_data;
	Parameter snpnames;
	Parameter traitnames;
	Parameter traits;
	Parameter snps;
	Parameter heritabilities;
	Parameter chi;
	Parameters();
	Parameters(int, char*[]);		//	Constructor from cmdline
	static string get_cmd_line(int,char*[]);
	bool defaultstate;//Check: are there parameters instead of infile's path

};

ostream &operator <<(ostream &, Parameters);

#endif /* PARAMETERS_H_ */
