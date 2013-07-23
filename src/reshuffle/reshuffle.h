/*
 * reshuffle.h
 *
 *  Created on: 01.03.2013
 *      Author: lima
 */

#ifndef RESHUFFLE_H_
#define RESHUFFLE_H_

typedef long long int64_t;

#include "Parameters.h"
#include "iout_file.h"

using namespace std;
typedef long long int64_t;
class Reshuffle{
public:
	iout_file * p_iout_file;
	Parameters * p_Parameters;
	int per_trait_per_snp;
	set<int> traits2write;
	set<int> snps2write;
	string create_filename(string, string);
	string create_filename(string);
	void write_datadims(ofstream&);
	void write_snpnames(ofstream&);
	void write_traitnames(ofstream&);
	void write_data(ifstream&,ofstream&);
	void write_data_chi(ifstream&,ofstream&);
	void write_slim_data(ifstream&,ofstream&);
	void write_herest(ifstream&,ofstream&);
	int64_t herest_startpos;
	int est_shift(int);
	int est_beta_shift(int);
	void run(Parameters);
	Reshuffle(iout_file &,Parameters &);
	void run();
};
#endif /* RESHUFFLE_H_ */
