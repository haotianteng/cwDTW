#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/proc.h"
#include "util/exception.h"
#include "kmer/kmer_index.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <iomanip>

using namespace std;
using namespace g::proc;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot);
void getRootName(string &in,string &out,char slash);
//---------------------- utility ------//over
//--------- pore_models: from genome to expected signal ----------//
bool Genomes2SignalSequence(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals, int scale, int FIVE_or_SIX);


//--------- pore_models: from genome to expected signal (RNA) ----------//
bool Genomes2SignalSequence_RNA(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals, int scale, int mv200_or_mv180);

//--------------- continuous wavelet transform (CWT) analysis -----------------//
/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void CWTAnalysis(const std::vector<double>& raw, std::vector<std::vector<double> >& output, 
	double scale0, double dscale, int npyr);

//----------------- boundary generation for constrained dynamic time warping (cDTW) -------------//
void BoundGeneration(std::vector<std::pair<int,int> >& cosali, 
	int neib, std::vector<std::pair<int,int> >& bound, int mode);

//====================== FastDTW ===================//
void FastDTW(std::vector<double>& in1, std::vector<double>& in2,
	std::vector<std::pair<int,int> >& alignment, int radius, int max_level, int mode,
	double* totaldiff = 0);


//====================== continous wavelet dynamic time warping (cwDTW) ========================//
void MultiLevel_WaveletDTW(std::vector<double>& in1, std::vector<double>& in2,
	std::vector<std::vector<double> >& sig1, std::vector<std::vector<double> >& sig2, 
	std::vector<std::pair<int,int> >& alignment, int radius, int test, int mode, 
	double* totaldiff = 0);

//----------- write alignment to file (nano) -------------//__2017.10.15__//(Sheng modified)
void WriteSequenceAlignment_nano(const char* output,int KMER, 
	const std::vector<double>& reference, const std::vector<double>& peer,
	const std::vector<int>& refer_orig, const std::vector<int>& peer_orig,
	vector<pair<int,int> >& alignment, int swap, std::string &refer_str);


//------------------ write alignment to file ----------------------//
void WriteSequenceAlignment(const char* output,
        const std::vector<double>& reference, const std::vector<double>& peer,
        vector<pair<int,int> >& alignment, int swap);