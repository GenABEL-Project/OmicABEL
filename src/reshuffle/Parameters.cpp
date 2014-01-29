/*
 * Parameters.cpp
 *
 *  Created on: 28.02.2013
 *      Author: Sodbo
 */
#include "Parameters.h"

using namespace std;

string sep = "--"; // Set separator for cmdline

//Default constructor
Parameter::Parameter() {
	name = "None";
	use = false;
	def_value = true;
}

//Overloading operator cout for class Parameter
ostream &operator <<(ostream &os, Parameter par) {
	os << "\nPARAMETR" << "\t" << par.name << "\t" << "USE" << "\t" << par.use<<"\t"<<"Use def_val "<<par.def_value;
	cout<<"\tNumbers set ";
	for (set<int>::iterator it= par.numbersset.begin();it!=par.numbersset.end();it++)
		os <<*it<<" ";
	cout<<"\tNames set ";
	for (set<string>::iterator it= par.namesset.begin();it!=par.namesset.end();it++)
		os <<*it<<" ";
	return os;
}

//Overloading operator cout for class Parameters
ostream &operator <<(ostream &os, Parameters par) {
	os << "\nOmicABEL iout_file is "<<par.iout_fname;
	os << "\nOmicABEL out_file is "<<par.out_fname;
	if (par.default_outfile){
		os << "\nUsing default out file name";
	}else
		os << "\nOut file name : "<<par.outfile;
	if(par.get_help)
		os << "\nPrinting help";
	if(par.get_info)
		os << "\nPrinting info about run";
	if(par.write_datadims)
		os << "\nWriting data dimension";
	if(par.snpnames.use)
		os << par.snpnames;
	if(par.traitnames.use)
		os << par.traitnames;
	if(par.traits.use)
		os << par.traits;
	if(par.snps.use)
		os << par.snps;
	if(par.herit.use)
		os << par.herit;
	if(par.chi.use)
		os << par.chi;
	if(par.write_slim_data)
			os << "\nWriting slim data";
	return os;
}

// default constructor
Parameters::Parameters(){
}

//Constructor-Parser from command line argument
Parameter::Parameter(char* cmd_val, string paramname) {
	name = paramname;
	use=true;
	if (cmd_val!=NULL) {
		def_value=false;
		string value = cmd_val;
		string val_d=value+",";
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
						value.erase(value.find("file="),5+str_tmp.size());//kill filename substring in value
					}else
						namesset.insert(str_tmp);
				}
				str_tmp="";
			}
		}
	} else{
		def_value = true; // default value for parameters is "all"
	}
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
	if (filename.find("--help")==0||filename.find("-h")==0)
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
				traitnames = Parameter(optarg,"traitnames");
				break;
			}
			case 'm':{
				snpnames = Parameter(optarg,"snpnames");
				break;
			}
			case 't':{
				traits = Parameter(optarg,"traits");
				break;
			}
			case 's':{
				snps = Parameter(optarg,"snps");
				break;
			}
			case 'e':{
				herit = Parameter(optarg,"herit");
				break;
			}
			case 'c':{
				chi = Parameter(optarg,"chi");
				if (chi.numbersset.size()==0){
					cout<<"\nChi value is not a number. Please, set correct Chi value";
					exit(1);
				}

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
		defaultstate=write_datadims+snpnames.use+traitnames.use+traits.use
					+snps.use+herit.use+chi.use+write_slim_data;
		param_coutner = write_datadims+snpnames.use+traitnames.use +
				(traits.use||snps.use||(chi.use&&!write_slim_data)) +
				(write_slim_data&&chi.use);
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
