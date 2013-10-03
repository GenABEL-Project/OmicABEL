/*
 * reshuffle.cpp
 *
 *  Created on: 01.03.2013
 *      Author: lima
 */

#include "reshuffle.h"

using namespace std;

#define PRECISION_DOUBLE 15//Precision of double

Reshuffle::Reshuffle(iout_file &iout,Parameters &Params){
	p_iout_file = &iout;
	p_Parameters = &Params;
	per_trait_per_snp = (*p_iout_file).header.p + (*p_iout_file).header.p * ((*p_iout_file).header.p + 1) / 2;
	herest_startpos = (long long)(p_iout_file->header.m * p_iout_file->header.t*
			(p_iout_file->header.p + p_iout_file->header.p * (p_iout_file->header.p + 1) / 2))
					* sizeof(double);
}

void Reshuffle::write_datadims(ofstream& txt_datadims){

	txt_datadims << "Number of traits\t" << (*p_iout_file).header.t << endl;
	txt_datadims << "Number of SNP\t" << (*p_iout_file).header.m << endl;
	txt_datadims << "Number of covariates\t" << ((*p_iout_file).header.p - 2);
	cout<<"End write data dimension\t"<< double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_snpnames(ofstream& txt_snpnames){

	if ((*p_Parameters).snpnames.def_value){
		for (unsigned int i=0;i<(*(*p_iout_file).labels.snp_names).size();i++)
			(*p_Parameters).snpnames.numbersset.insert(i);
		cout<<"\nWriting all SNP's names";
	}
	for(set<int>::iterator it= (*p_Parameters).snpnames.numbersset.begin();it!=(*p_Parameters).snpnames.numbersset.end();++it)
		txt_snpnames << "\nSNP #"<<(*it+1)<<"\t"<<(*(*p_iout_file).labels.snp_names)[*it];
	cout<<"\nEnd write SNP's names\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";
}

void Reshuffle::write_traitnames(ofstream& txt_traitnames){

	if ((*p_Parameters).traitnames.def_value){
		for (unsigned int i=0;i<(*(*p_iout_file).labels.trait_names).size();i++)
			(*p_Parameters).traitnames.numbersset.insert(i);
		cout<<"\nWriting all trait's names";
	}
	for(std::set<int>::iterator it= (*p_Parameters).traitnames.numbersset.begin();it!=p_Parameters->traitnames.numbersset.end();++it)
		txt_traitnames<<"TRAIT #"<<(*it+1)<<"\t"<<(*(*p_iout_file).labels.trait_names)[*it];
	cout<<"\nEnd write trait's names\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_data(ifstream& out_file,ofstream& data){
	out_file.seekg(0, ios_base::beg);
	cout << "\nStart write data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";
	data.precision(PRECISION_DOUBLE);
	data<<"SNP\t";
	data<<	"Trait\t";
	for (unsigned int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
		data << (*(*p_iout_file).labels.beta)[beta] << "\t";
	for (unsigned int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
		data << (*(*p_iout_file).labels.se)[se] << "\t";
	for (unsigned int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
		data << (*(*p_iout_file).labels.cov)[cov] << "\t";
	data << endl;
	double* buf = new double[per_trait_per_snp];
	char s[30];
	ostringstream ostr;
	for (set<int>::iterator trait= (*p_Parameters).traits.numbersset.begin();trait!=(*p_Parameters).traits.numbersset.end();trait++) {
		long long oldPos=0,pos = 0;
		for (set<int>::iterator snp= (*p_Parameters).snps.numbersset.begin();snp!=(*p_Parameters).snps.numbersset.end();snp++) {
			ostr << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
			ostr << (*(*p_iout_file).labels.trait_names)[*trait]<<"\t";
			pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			for (int i = 0; i < per_trait_per_snp; i++) {
				sprintf(s, "%.15g", buf[i]);
				ostr << s << "\t";
			}
			ostr << endl;
		}
		data << ostr.str();
		ostr.str("");
		ostr.clear();
	}
	delete buf;
	cout << "\nFinish write data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";
}

void Reshuffle::write_data_chi(ifstream& out_file,ofstream& txt_chi){
	out_file.seekg(0, ios_base::beg);
	double chi = 0;
	set<int>::iterator chi_val = (*p_Parameters).chi.numbersset.begin();
	double CheckChi = *chi_val+1;
	cout << "\nStart write chi data=" << double(clock()) / CLOCKS_PER_SEC <<" sec";
	txt_chi.precision(PRECISION_DOUBLE);
	txt_chi << "SNP\t";
	txt_chi << "Trait\t";
	for (unsigned int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
		txt_chi << (*(*p_iout_file).labels.beta)[beta] << "\t";
	for (unsigned int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
		txt_chi << (*(*p_iout_file).labels.se)[se] << "\t";
	for (unsigned int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
		txt_chi << (*(*p_iout_file).labels.cov)[cov] << "\t";
	txt_chi << "Chi2" << endl;
	double* buf = new double[per_trait_per_snp];
	char s[30];
	ostringstream ostr;
	for (set<int>::iterator trait= (*p_Parameters).traits.numbersset.begin();trait!=(*p_Parameters).traits.numbersset.end();trait++) {
		long long oldPos=0,pos = 0;
		for (set<int>::iterator snp= (*p_Parameters).snps.numbersset.begin();snp!=(*p_Parameters).snps.numbersset.end();snp++) {
			pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			chi=pow((buf[(*(*p_iout_file).labels.beta).size()-1]/buf[(*(*p_iout_file).labels.beta).size()+(*(*p_iout_file).labels.se).size()-1]),2);
			if(chi>CheckChi){
				ostr << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
				ostr << (*(*p_iout_file).labels.trait_names)[*trait]<<"\t";
				for (int i = 0; i < per_trait_per_snp; i++) {
					sprintf(s, "%.15g", buf[i]);
					ostr << s << "\t";
				}
				ostr << chi << endl;
			}
		}
		txt_chi << ostr.str();
		ostr.str("");
		ostr.clear();
	}
	delete buf;
	cout << "\nFinish write chi data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";
}

void Reshuffle::write_slim_data(ifstream& out_file, ofstream& txt_slim){
	cout << "\nStart write slim data=" << double(clock()) / CLOCKS_PER_SEC <<" sec";
	out_file.seekg(0, ios_base::beg);
	set<int> goodtraits;
	set<int> goodsnps;
	double chi = 0;
	if((*p_Parameters).chi.def_value||(!(*p_Parameters).chi.use)){
		cout << "\nERROR: " << "Chi value doesn't set";
		cout << "\nPlease, set Chi value to get slim data";
		exit(1);
	}
	set<int>::iterator chi_val = (*p_Parameters).chi.numbersset.begin();
	double CheckChi = *chi_val+1;
	double* buf = new double[per_trait_per_snp];
	for (set<int>::iterator trait= (*p_Parameters).traits.numbersset.begin();trait!=(*p_Parameters).traits.numbersset.end();trait++) {
		long long oldPos=0,pos = 0;
		for (set<int>::iterator snp= (*p_Parameters).snps.numbersset.begin();snp!=(*p_Parameters).snps.numbersset.end();snp++) {
			pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			//cout << oldPos << "-" << pos << endl;
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			chi=pow((buf[(*(*p_iout_file).labels.beta).size()-1]/buf[(*(*p_iout_file).labels.beta).size()+(*(*p_iout_file).labels.se).size()-1]),2);
			if(chi>CheckChi){
				goodtraits.insert(*trait);
				goodsnps.insert(*snp);
			}
		}
	}
	ostringstream ostr;
	out_file.seekg(0, ios_base::beg);
	ostr << "SNP\t";
	ostr << "Trait\t";
	for (unsigned int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
		ostr << (*(*p_iout_file).labels.beta)[beta] << "\t";
	for (unsigned int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
		ostr << (*(*p_iout_file).labels.se)[se] << "\t";
	for (unsigned int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
		ostr << (*(*p_iout_file).labels.cov)[cov] << "\t";
	ostr << "Chi2" << endl;

	char s[30];
	for (set<int>::iterator trait= goodtraits.begin();trait!=goodtraits.end();trait++) {
		long long oldPos=0,pos = 0;
		for (set<int>::iterator snp= goodsnps.begin();snp!=goodsnps.end();snp++) {
			ostr << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
			ostr << (*(*p_iout_file).labels.trait_names)[*trait]<<"\t";
			pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			chi=pow((buf[(*(*p_iout_file).labels.beta).size()-1]/buf[(*(*p_iout_file).labels.beta).size()+(*(*p_iout_file).labels.se).size()-1]),2);
			for (int i = 0; i < per_trait_per_snp; i++) {
				sprintf(s, "%.15g", buf[i]);
				ostr << s << "\t";
			}
			ostr << chi << endl;
		}
		txt_slim << ostr.str();
		ostr.str("");
		ostr.clear();
	}
	delete buf;
	cout <<"\nEnd write slim data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";
}

int Reshuffle::est_shift(int counter){
	int shift = ((p_iout_file->header.m * p_iout_file->header.t* (p_iout_file->header.p +
			p_iout_file->header.p * (p_iout_file->header.p + 1) / 2))
			+counter*p_iout_file->header.t)* sizeof(double);
	return shift;
}

int Reshuffle::est_beta_shift(int counter){
	int shift = est_shift(3)+counter*(p_iout_file->header.p-1)*sizeof(double);
	return shift;
}

void Reshuffle::write_herest(ifstream& out_file, ofstream& herest){
	cout << "\nStart write heritabilities and estimates=" << double(clock()) / CLOCKS_PER_SEC <<" sec";
	ofstream txt_est("estimates.txt");
	out_file.seekg(herest_startpos, ios_base::beg);
	if (p_Parameters->herit.def_value)
		for(unsigned int i=0;i<(*(p_iout_file->labels.trait_names)).size();i++)
			p_Parameters->herit.numbersset.insert(i);
	txt_est.precision(PRECISION_DOUBLE);
	txt_est<<"\t";
	for (set<int>::iterator trait= p_Parameters->herit.numbersset.begin();trait!=p_Parameters->herit.numbersset.end();trait++)
		txt_est << (*(p_iout_file->labels.trait_names))[*trait] << "\t";
	txt_est << endl;
	list<string> est_names;
	est_names.insert(est_names.end(), "h2");
	est_names.insert(est_names.end(), "var");
	est_names.insert(est_names.end(), "resVar");
	double tmp_number = 0;
	int counter=0;
	for (list<string>::iterator name = est_names.begin();name != est_names.end(); ++name) {
		txt_est << *name << "\t";
		for (std::set<int>::iterator trait= p_Parameters->herit.numbersset.begin();trait!=p_Parameters->herit.numbersset.end();++trait) {
			out_file.seekg(*trait*sizeof(double),ios_base::cur);
			out_file.read((char *) &tmp_number, sizeof(double));
			txt_est << tmp_number << "\t";
			out_file.seekg(est_shift(counter), ios_base::beg);
		}
		counter++;
		txt_est << "\n";
		out_file.seekg(est_shift(counter), ios_base::beg);
	}
	out_file.seekg(est_shift(3), ios_base::beg);
	counter=0;
	for (unsigned int beta=0;beta<(*(p_iout_file->labels.beta)).size();beta++) {
		beta++;
		if (beta != (*(p_iout_file->labels.beta)).size()) {
			beta--;
			txt_est << (*(p_iout_file->labels).beta)[beta] << "\t";
			for (std::set<int>::iterator trait= p_Parameters->herit.numbersset.begin();trait!=p_Parameters->herit.numbersset.end();++trait) {
				out_file.seekg(*trait*sizeof(double),ios_base::cur);
				out_file.read((char *) &tmp_number, sizeof(double));
				txt_est << tmp_number << "\t";
				out_file.seekg(est_beta_shift(counter), ios_base::beg);
			}
			counter++;
			out_file.seekg(est_beta_shift(counter), ios_base::beg);
			txt_est << "\n";
			beta++;
		}
	}
	cout << "\nEnd write heritabilities and estimates\t" << double(clock()) / CLOCKS_PER_SEC <<" sec";
}

void Reshuffle::run(){
	if((*p_Parameters).write_datadims){
		ofstream datadims((*p_Parameters).outfile.c_str());
		write_datadims(datadims);
	}
	if((*p_Parameters).snpnames.use){
		ofstream snpnames((*p_Parameters).outfile.c_str());
		write_snpnames(snpnames);
	}

	if((*p_Parameters).traitnames.use){
		ofstream traitnames((*p_Parameters).outfile.c_str());
		write_traitnames(traitnames);
	}

	//Open *.out to read data
	ifstream out_file((*p_Parameters).out_fname.c_str(), ios::binary | ios::in);
	if (!out_file) {
		cout << "\nError open " << (*p_Parameters).out_fname << " file";
		cout << "\nMaybe, file doesn't exist";
		exit(1);
	}

	//If any of parameters traits||snps||chi use, this block fill traits.numbersset and snps.numbersset
	//(if their values are default all)
	if((*p_Parameters).traits.use||(*p_Parameters).snps.use||(*p_Parameters).chi.use||!(*p_Parameters).defaultstate){

			if((*p_Parameters).traits.def_value||(!(*p_Parameters).traits.use)){
				for(unsigned int i=0;i<(*(*p_iout_file).labels.trait_names).size();i++)
				(*p_Parameters).traits.numbersset.insert(i);
			cout<<"\nWriting data for all traits";
		}

		if((*p_Parameters).snps.def_value||(!(*p_Parameters).snps.use)){
			for(unsigned int i=0;i<(*(*p_iout_file).labels.snp_names).size();i++)
				(*p_Parameters).snps.numbersset.insert(i);
			cout<<"\nWriting data for all SNPs";
		}
	}

	if((((*p_Parameters).traits.use||(*p_Parameters).snps.use)&&!(*p_Parameters).chi.use)||!(*p_Parameters).defaultstate){
		ofstream data((*p_Parameters).outfile.c_str());
		write_data(out_file,data);
	}

	if((*p_Parameters).chi.use&&!(*p_Parameters).write_slim_data){
		ofstream chi_data((*p_Parameters).outfile.c_str());
		write_data_chi(out_file,chi_data);

	}

	if((*p_Parameters).write_slim_data){
		ofstream dataslim((*p_Parameters).outfile.c_str());
		write_slim_data(out_file,dataslim);

	}

	if((*p_Parameters).herit.use){
		ofstream herest((*p_Parameters).outfile.c_str());
		write_herest(out_file,herest);
	}
}
