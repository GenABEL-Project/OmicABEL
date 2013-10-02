/*
 * Parameters.cpp
 *
 *  Created on: 28.02.2013
 *      Author: Sodbo
 */
#include "Parameters.h"

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

//Overloading operator cout for class Parameter
ostream &operator <<(ostream &os, Parameter par) {
	os << "PARAMETR" << "\t" << par.name << "\t" << "USE" << "\t" << par.use
			<< "\t" << "VALUE" << "\t" << par.value;
	cout<<"\tNumbers set ";
	for (set<int>::iterator it= par.numbersset.begin();it!=par.numbersset.end();it++)
		os <<*it<<" ";
	cout<<"\tNames set ";
	for (set<string>::iterator it= par.namesset.begin();it!=par.namesset.end();it++)
		os <<*it<<" ";
	//os<<"\toutfile="<<par.outfile;
	os<<endl;
	return os;
}

//Overloading operator cout for class Parameters
ostream &operator <<(ostream &os, Parameters par) {
	os << "OmicABEL iout_file is "<<par.iout_fname<<endl;
	os << "OmicABEL out_file is "<<par.out_fname<<endl;
	if (par.default_outfile){
		os << "Using default out file name"<<endl;
	}else
		os << "Out file name : "<<par.outfile<<endl;

	os << "Print help : " << par.get_help<<endl;
	os << "Print info : " << par.get_info<<endl;
	os << "Write data dimension : " << par.write_datadims<<endl;
	os << par.snpnames;
	os << par.traitnames;
	os << par.traits;
	os << par.snps;
	os << par.heritabilities;
	os << par.chi;
	os << "Write slim data : " << par.write_slim_data;
	return os;
}

// default constructor
Parameters::Parameters(){
}

//Constructor-Parser from command line
Parameter::Parameter(string cmdline, string paramname, string ofile) {

	name = paramname;
	outfile=ofile;
	int parpos = cmdline.find("--"+name)+2; //	position of substring with param's name
	if (parpos != string::npos+2) { //	check  [FOUND OR NOT]
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
							numbersset.insert(i);
						}
					}else if(atoi(str_tmp.c_str())!=0){
						numbersset.insert(atoi(str_tmp.c_str())-1);
					}else{
						if(str_tmp.find("file=")!=string::npos){//Fill outfile name
							str_tmp=str_tmp.erase(0,5);
							outfile=str_tmp;
							value.erase(value.find("file="),5+str_tmp.size());//kill filename substring in value
						}else
							namesset.insert(str_tmp);
					}
					str_tmp="";
				}
			}
		} else
			value = "all"; // default value for parameters is "all"
		if(value.size()==0)
			value="all";

	} else { // If parameter aren't in cmdline
		use = false;
		value = "None";
	}
}

string Parameter::delfromcmdline(string cmdline){
	if(use){
		int val=value.size();
		if(value=="None"||value=="all")
			val=-1;
		cmdline.replace(cmdline.find(name)-2,name.size()+val+3,"");
	}
	return cmdline;
}

//	Constructor, which gets info about all parameters from cmdline
Parameters::Parameters(int argc,char* argv[]) {

	get_help = false;
	get_info = false;
	write_datadims = false;
	write_slim_data = false;
	default_outfile = true;

	//get input file names and delete it from argv[]
	string filename= argv[1];
	if (filename.find("--")==0||filename.find("-h")==0)
		get_help=true;
	iout_fname = filename + ".iout";
	out_fname = filename + ".out";
	argv[1] = " ";

	const char* short_options = "hidn::m::t::s::e::c:l:o";

	const struct option long_options[] = {
			{"help",no_argument,NULL,'h'},
			{"info",no_argument,NULL,'i'},
			{"datadims",no_argument,NULL,'d'},
			{"traitnames",optional_argument,NULL,'n'},
			{"snpnames",optional_argument,NULL,'m'},
			{"traits",optional_argument,NULL,'t'},
			{"snps",optional_argument,NULL,'s'},
			{"herit",optional_argument,NULL,'e'},
			{"chi",required_argument,NULL,'c'},
			{"dataslim",no_argument,NULL,'l'},
			{"outfile",required_argument,NULL,'o'},
			{NULL,0,NULL,0}
	};

	int rez = 0;
	int option_index = 0;
		while((rez=getopt_long(argc,argv,short_options,long_options,&option_index))!=-1){
			switch (rez){
			case 'h' :{
				get_help=true;
				break;
			}
			case 'i':{
				get_info=true;
				break;
			}
			case 'd':{
				write_datadims=true;
				break;
			}
			case 'n':{
				write_datadims=true;
				break;
			}
			case 'm':{
				write_datadims=true;
				break;
			}
			case 't':{
				write_datadims=true;
				break;
			}
			case 's':{
				write_datadims=true;
				break;
			}
			case 'e':{
				write_datadims=true;
				break;
			}
			case 'c':{
				write_datadims=true;
				break;
			}
			case 'l':{
				write_slim_data=true;
				break;
			}
			case 'o':{
				default_outfile = false;
				outfile=optarg;
				break;
			}
			}

		}
	//
	static string cmdline=get_cmd_line(argc,argv);
/*	snpnames = Parameter(cmdline, "snpnames","snpnames.txt");
	cmdline=snpnames.delfromcmdline(cmdline);
	traitnames = Parameter(cmdline, "traitnames","traitnames.txt");
	cmdline=traitnames.delfromcmdline(cmdline);
	traits = Parameter(cmdline, "traits","data.txt");
	cmdline=traits.delfromcmdline(cmdline);
	snps = Parameter(cmdline, "snps","data.txt");
	cmdline=snps.delfromcmdline(cmdline);
	heritabilities = Parameter(cmdline, "heritabilities","estimates.txt");
	cmdline=heritabilities.delfromcmdline(cmdline);
	chi = Parameter(cmdline, "chi","chi_data.txt");
	cmdline=chi.delfromcmdline(cmdline);
	defaultstate=write_datadims+snpnames.use+traitnames.use+traits.use
			+snps.use+heritabilities.use+chi.use+write_slim_data+run_test;
	if(traits.outfile!="data.txt"&&snps.outfile!="data.txt"){
		cout<<"You've set outfile name in <<traits>> and <<snps>> parameters"<<endl;
		cout<<"Please, set outfile name for data ones"<<endl;
		exit(1);
	}*/

}

void Parameter::setbynames(vector<string> names){
	int find=0;
	int before=0;
	int after=0;
	string regexp="";
	if(name=="snps"){
		for (set<string>::iterator nameit=namesset.begin();nameit!=namesset.end();++nameit){
			string tmp=*nameit;
			if(tmp.find("before")!=string::npos){
				tmp.erase(0,7);
				before = atoi(tmp.c_str());
			}
			if(tmp.find("after")!=string::npos){
				tmp.erase(0,6);
				after = atoi(tmp.c_str());
			}
			if(((*nameit).find("before")==string::npos)&&((*nameit).find("after")==string::npos)){
				for(int i=0;i<names.size();i++){
					if(names[i]==*nameit){
						for(int k=(i-before);k<(i+after+1);k++){
							numbersset.insert(k);
						}
						find=1;
						break;
					}
				}
				if(find==1)
					continue;
			}
		}
	}
	if(name=="traits"||name=="heritabilities"){
		for (set<string>::iterator nameit=namesset.begin();nameit!=namesset.end();++nameit){
			string tmp=*nameit;
			if(tmp.find("regexp")!=string::npos){
				tmp.erase(0,7);
				regexp = tmp;
				cout<<"REGEXP="<<regexp<<endl;
			}
		}
			for (set<string>::iterator nameit=namesset.begin();nameit!=namesset.end();++nameit){
				if((*nameit).find("regexp")==string::npos){
					for(int i=0;i<names.size();i++){
						if(names[i]==*nameit){
							numbersset.insert(i);
							find=1;
							break;
						}
					}
					if(find==1)
						continue;
				}else{
					for(int i=0;i<names.size();i++){
						if((names[i]).find(regexp)==0){
							numbersset.insert(i);
							find=1;
						}
					}
				}
			}
		}
}
