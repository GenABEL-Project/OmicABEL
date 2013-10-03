/*
 * main.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Sodbo Sharapov sharapovsodbo@gmail.com
 */
#include "iout_file.h"
#include "Parameters.h"
#include "reshuffle.h"


using namespace std;

void print_help(){
	cout <<endl;
	cout << "	Available commands"<<endl<<endl;
	cout << " --datadims | to get data's dimension"<<endl<<endl;
	cout << " --snpnames=<indexes> | to get names of SNP by indexes"<<endl<<endl;
	cout << " --traitnames=<indexes> | to get names of trait by indexes"<<endl<<endl;
	cout << " --traits=<indexes OR/AND names OR/AND regexp> | to get data"<<endl;
	cout << "   by trait's indexes OR/AND names OR/AND regexp"<<endl<<endl;
	cout << " --snp=<indexes OR/AND names OR/AND focus> | to get data"<<endl;
	cout << "   by SNP's indexes OR/AND names"<<endl<<endl;
	cout << " --heritabilities=<indexes OR/AND names OR/AND regexp>"<<endl;
	cout << "   to get estimates of trait's heritability, sigma, res_sigma and betas"<<endl<<endl;
	cout << " --chi=<number> | to get data by snp's, which chi2>number"<<endl<<endl;
	cout << " --dataslim | to get slim data [You should set --chi=<number>]"<<endl<<endl;
	exit(1);
}

int main(int argc, char* argv[]) {
	if (argc==1){
		print_help();
	}
	Parameters Params(argc, argv);
	if((Params.get_help||(Params.iout_fname==".iout"&&Params.out_fname==".out"))&&(!Params.run_test)){
		print_help();
	}
	if(Params.get_info)
		cout << Params;
	if(Params.param_coutner>1){
		cout<<"\nToo many parameters. Please, reduce number of parameters";
		exit(1);
	}

	cout<<"\nStart reshuffeling";
	iout_file iout_F(Params);
	cout << "\nFinish iout_file read\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";

	if(Params.traits.use)
		Params.traits.setbynames(*(iout_F.labels.trait_names));
	if(Params.snps.use)
		Params.snps.setbynames(*(iout_F.labels.snp_names));
	if(Params.herit.use)
		Params.herit.setbynames(*(iout_F.labels.trait_names));
	Reshuffle reshh(iout_F,Params);
	reshh.run();
	cout << "\nFinish reshuffling " << double(clock()) / CLOCKS_PER_SEC <<" sec";
	return (0);
}
