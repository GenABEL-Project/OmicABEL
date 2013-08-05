/*
 * main.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Sodbo
 */
#include <iostream>
#include <ctime>
#include "iout_file.h"
#include "Parameters.h"
#include "reshuffle.h"
#include <iterator>

using namespace std;

int main(int argc, char* argv[]) {
	//cout<<"Start reshuffeling"<<endl;
	Parameters Params(argc, argv);
	if((Params.help.use||(Params.iout_fname==".iout"&&Params.out_fname==".out"))&&(!Params.test.use)){
		cout<<"Available commands"<<endl<<endl;
		cout<<" --datadims | to get data's dimension"<<endl<<endl;
		cout<<" --snpnames=<indexes> | to get names of SNP by indexes"<<endl<<endl;
		cout<<" --traitnames=<indexes> | to get names of trait by indexes"<<endl<<endl;
		cout<<" --traits=<indexes OR/AND names OR/AND regexp> | to get data"<<endl;
		cout<<"   by trait's indexes OR/AND names OR/AND regexp"<<endl<<endl;
		cout<<" --snp=<indexes OR/AND names OR/AND focus> | to get data"<<endl;
		cout<<"   by SNP's indexes OR/AND names"<<endl<<endl;
		cout<<" --heritabilities=<indexes OR/AND names OR/AND regexp>"<<endl;
		cout<<"   to get estimates of trait's heritability, sigma, res_sigma and betas"<<endl<<endl;
		cout<<" --chi=<number> | to get data by snp's, which chi2>number"<<endl<<endl;
		cout<<" --dataslim | to get slim data [You should set --chi=<number>]"<<endl<<endl;
		exit(1);
	}

	if(Params.info.use)
		cout << Params;
	iout_file iout_F(Params);
	cout << "Finish iout_file read\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	if(Params.info.use){
		cout<<iout_F.header;
		cout<<iout_F.labels;
	}
	if(Params.traits.use)
		Params.traits.setbynames(*(iout_F.labels.trait_names));
	if(Params.snps.use)
		Params.snps.setbynames(*(iout_F.labels.snp_names));
	if(Params.heritabilities.use)
		Params.heritabilities.setbynames(*(iout_F.labels.trait_names));
	Reshuffle reshh(iout_F,Params);
	reshh.run();
	cout << "Finish reshuffling " << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	return (0);
}
