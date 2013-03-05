/*
 * Parameters.cpp
 *
 *  Created on: 28.02.2013
 *      Author: lima
 */
#include "Parameters.h"
#include <stdlib.h>
#include <string>
using namespace std;

string sep = "--"; // Set separator for cmdline

//	Convert cmdline to string
string Parameters::get_cmd_line(int argc, char* argv[]) {
	string cmd_line = "";
	for (int i = 1; i < argc; i++) {
		cmd_line += argv[i];
	}
	cmd_line += sep; // Need for cmdline parsing
	return cmd_line;
}
//Default constructor
Parameter::Parameter() {
	name = "None";
	use = false;
	value = "all";
}

ostream &operator <<(ostream &os, Parameter par) {
	os << "PARAMETR" << "\t" << par.name << "\t" << "USE" << "\t" << par.use
			<< "\t" << "VALUE" << "\t" << par.value;
	cout<<"\tValue set ";
	for (set<int>::iterator it= par.valueset.begin();it!=par.valueset.end();it++)
		os <<*it<<" ";
	os<<endl;
	return os;
}

ostream &operator <<(ostream &os, Parameters par) {
	os << "IOUT file is "<<par.iout_fname<<endl;
	os << "OUT file is "<<par.out_fname<<endl;
	os << par.datadims;
	os << par.snpnames;
	os << par.traitnames;
	os << par.traits;
	os << par.snps;
	os << par.heritabilities;
	os << par.chi;
	os << par.dataslim;
	return os;
}

//Constructor-parser from command line
Parameter::Parameter(string cmdline, string paramname) {

	name = paramname;
	int parpos = cmdline.find(name); //	position of substring with param's name
	if (parpos != string::npos) { //	check  [FOUND OR NOT]
		string val = "";
		unsigned int iter = parpos + name.size() + 1; //	iterator, which run from "=" to sep
		unsigned int seppos = cmdline.find(sep, parpos); // 	separator position
		while (iter < seppos) {
			val += cmdline.at(iter); // build value of parameter
			iter++;
		}
		use = true;
		if (val.size() != 0) {
			value = val;
			string val_d=val+",";
			string str_tmp="";
			for(unsigned int i=0;i<val_d.size();i++){
				if(val_d[i]!=','){
				str_tmp+=val_d[i];
				}else {
					string str2="";
					int postire=str_tmp.find("-");
					if(postire!=string::npos){
						string start="";
						string end="";
						for(int i=0;i<postire;i++){
							start+=str_tmp[i];
						}
						for(int i=(postire+1);i<str_tmp.size();i++){
							end+=str_tmp[i];
						}
						for(int i=atoi(start.c_str())-1;i<atoi(end.c_str());i++){
							valueset.insert(i);
						}
					}else valueset.insert(atoi(str_tmp.c_str())-1);
					str_tmp="";
				}
			}
		} else
			value = "all"; // default value for parameters is "all"
	} else { // If parameter aren't in cmdline
		use = false;
		value = "None";
	}
}

//	Constructor, which gets info about all parameters from cmdline
Parameters::Parameters(int argc,char* argv[]) {
	string cmdline=get_cmd_line(argc,argv);
	unsigned int seppos = cmdline.find(sep); // 	first separator's position
	string filename = cmdline.substr(0, seppos);
	iout_fname = filename + ".iout";
	out_fname = filename + ".out";
	datadims = Parameter(cmdline, "datadims");
	snpnames = Parameter(cmdline, "snpnames");
	traitnames = Parameter(cmdline, "traitnames");
	traits = Parameter(cmdline, "traits");
	snps = Parameter(cmdline, "snps");
	heritabilities = Parameter(cmdline, "heritabilities");
	chi = Parameter(cmdline, "chi");
	dataslim= Parameter(cmdline, "dataslim");
}
