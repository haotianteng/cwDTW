#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <proc/io.h>
#include "util/exception.h"
#include "5mer/5mer_index.h"

bool Genomes2SignalSequence(const std::vector<char>& genomes, std::vector<int>& signals, int scale)
{
	size_t bound = genomes.size()-5;//genomes.size()%5;
	for(size_t i = 0; i < bound; i++){
		int idx = g::Mer2Signal::FiveMer2Index(genomes[i], genomes[i+1], genomes[i+2], genomes[i+3], genomes[i+4]);
		double sigval = 3.8*g::Mer2Signal::AvgSignalAt(idx);
		
		for(int c = scale; c--;){
			signals.push_back((int)sigval);
		}
	}

//	for(size_t i = bound; i < genomes.size(); i++){
//		for(int c = scale; c--;){
//			signals.push_back(100);
//		}
//	}
}

int main(int argc, char **argv)
{
	std::string input="";
	std::string output="";

	struct options opts;
	opts.scale=1;
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	input=opts.input;
	output=opts.output;
	int scale=opts.scale;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		return -1;
	}

	std::vector<char> genomes;
	std::vector<int> signals;
	g::io::ReadATCG(opts.input, genomes);
	Genomes2SignalSequence(genomes, signals, scale);
	g::io::WriteSignalSequence_int(opts.output, signals);

	return 0;
}

