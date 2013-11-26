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
	cout << "\n	Available commands";
	cout << "\n\n --help | to get this help message";
	cout << "\n\n --datadims | to get data's dimension";
	cout << "\n\n --snpnames=<indexes> | to get names of SNP by indexes";
	cout << "\n\n --traitnames=<indexes> | to get names of trait by indexes";
	cout << "\n\n --traits=<indexes OR/AND names OR/AND regexp> | to get data";
	cout << "\n   by trait's indexes OR/AND names OR/AND regexp";
	cout << "\n\n --snp=<indexes OR/AND names OR/AND focus> | to get data";
	cout << "\n   by SNP's indexes OR/AND names";
	cout << "\n\n --herit=<indexes OR/AND names OR/AND regexp>";
	cout << "\n   to get estimates of trait's heritability, sigma, res_sigma and betas";
	cout << "\n\n --chi | to get data with chi-square column";
	cout << "\n\n --chi=<number> | to get data by snp's, which chi2>number";
	cout << "\n\n --dataslim | to get slim data with threshold over X";
	cout << "\n   [You should set --chi=<X>]";
	cout << "\n\n --info | to get some information about program's run"<<endl<<endl;
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
