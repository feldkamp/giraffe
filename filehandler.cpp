//
//  main.cpp
//  filehandler
//
//  Created by Feldkamp on 11/7/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <iomanip>

#include <string>
using std::string;

#include <sstream>

#include <vector>
using std::vector;

#include <fstream>

#include <boost/program_options.hpp>		//parser for options (commandline and ini-file)
namespace po = boost::program_options;

#include "arrayclasses.h"
#include "arraydataIO.h"




//-----------------------------------------------------------------------------------------------
string defaultDir(){
	string base = "";
	char *home = getenv( "HOME" );
	if (home){
		base = home;
		base += "/Desktop";
	}
	return base;
}


//-----------------------------------------------------------------------------------------------
int main (int argc, char * const argv[]) {

	arraydataIO *io = new arraydataIO;
	std::ofstream fout("filehandler.log");
	
	string mode = "";	
	vector<string> in_fn;
	string list_fn = "";
	string out_fn = "";
	double p1 = 0.;
	double p2 = 0.;

	po::positional_options_description p;
	p.add("mode", -1);

	po::options_description desc("Usage: filehandler <mode> [options] \nAllowed options:");
	desc.add_options()		
		("mode,m", 		po::value<string>(&mode),								"mode of operation,"
									"\nchoose a keyword:"
									"\n   add           (2 files)"
									"\n   sub           (2 files)"
									"\n   mult          (2 files)"
									"\n   div           (2 files)"
									"\n   mask_thresh   (1 file) --p1 --p2"
									"\n       p1 = lower threshold, p2 = upper threshold"
									"\n   mask_invert   (1 file)"
									"\n   mask_apply    (2 files)"
									"\n   stats         (x files)"
									"\n   avg           (x files)"
									"\n   IO            (x files)"
									"\n      output file type is determined by extension"
									)
		("param_one", 	po::value<double>(&p1), 								"generic parameter 1 (use depends on chosen mode)")
		("param_two", 	po::value<double>(&p2), 								"generic parameter 2 (use depends on chosen mode)")
		("output,o",	po::value<string>(&out_fn), 							"output filename")	
		("input,i", 	po::value< vector<string> >(&in_fn)->multitoken(),		"input filename")
		("list,l", 		po::value<string>(&list_fn),							"input filelist")
		("help,h",																"produce help message")	
	;
	
	po::variables_map vm;
	try{
		//po::parsed_options parsed = po::parse_command_line(argc, argv, desc);
		po::parsed_options parsed = po::command_line_parser(argc, argv).
				options(desc).positional(p).allow_unregistered().run();
		po::store(parsed, vm);
		po::notify(vm);
	}catch( const boost::program_options::error &e ){
		cerr << "Error in command line arguments: " << e.what() << "." << endl;
		cerr << desc << endl;
		return 1;
	}	
	//cout << ">parser done<" << endl;
	
	if (vm.count("help")) {
		cout << desc << endl;
		return 0;
	}
	if (vm.count("mode")){
		cout << "--> using mode " << mode << endl;
	}else{
		cerr << "Error. No mode of operation given" << endl;
		cerr << desc << endl;
		return 1;
	}
	if (vm.count("input")){
		for (int i = 0; i < in_fn.size(); i++){
			cout << "--> using input file " << in_fn.at(i) << endl;
		}
	}
	if (vm.count("list")){
		cout << "--> using input file list " << list_fn << endl;
		std::ifstream fin;
		fin.open(list_fn.c_str());					
		if( !fin.fail() ){
			string line;
			while (fin.good()) {
				getline(fin, line);
				if (line != ""){
					cout << line << endl;
					in_fn.push_back(line);
				}
			}
		}else{
			cerr << "Error. Could not open file list '" << list_fn << "', aborting..." << endl;
			return 2;
		}
		fin.close();
	}
	if (vm.count("output")){
		cout << "--> using output file " << out_fn << endl;
	}	
	
	if (out_fn == ""){
		out_fn = defaultDir() + "/out.h5";
	}


	unsigned int num_files = (unsigned int)in_fn.size();
	
	array2D *first = 0;
	array2D *second = 0;
	
	if (num_files == 0){
		cerr << "No input file given. " << endl;
		return 2;
	}
	

	//========================================================================== one-file operations
	if (mode == "mask_thresh" || mode == "mask_invert"){ 
		if (num_files > 0){
			int fail = io->readFromFile( in_fn.at(0), first );
			if (fail){
				cerr << "Error reading file" << endl;
				return 3;
			}
		}else{
			cerr << "Not enough input files for this mode" << endl;
			return 3;
		}

		if (mode == "mask_thresh"){
			cout << "THRESHOLDING" << endl;
			double lower = p1;
			double upper = p2;
			int keep, use_high, use_low = 0;
			array2D *mask_after_thresholding = new array2D( first->dim1(), first->dim2() );
			for (int i=0; i<first->size(); i++){
				if (first->get_atIndex(i) > lower){
					if (first->get_atIndex(i) <= upper){
						//keep original value
						mask_after_thresholding->set_atIndex(i, 1);
						keep++;
					} else {
						//use upper value
						first->set_atIndex(i, upper);
						mask_after_thresholding->set_atIndex(i, 0);
						use_high++;
					}
				} else {
					//use lower value
					first->set_atIndex(i, lower);
					mask_after_thresholding->set_atIndex(i, 0);
					use_low++;
				}
			}
			string mask_fn = "mask_after_thresholding.h5";
			io->writeToHDF5( mask_fn, mask_after_thresholding );
			delete mask_after_thresholding;
			cout << first->size() << " values. " << keep << " kept, " 
				<< use_high << " replaced by '" << upper << "' (high limit), "
				<< use_low << " replaced by '" << lower << "' (low limit)" << endl;
		}
		
		if (mode == "mask_invert"){
			cout << "INVERSION OF MASK (0-->1, 1-->0)" << endl;
			double checkval = 0.1;
			for (int i=0; i<first->size(); i++){
				if (first->get_atIndex(i) > checkval){
					first->set_atIndex(i, 0);
				} else {
					first->set_atIndex(i, 1);
				}
			}
		}		
		
		//write output
		io->writeToFile( out_fn, first );
	}//end of one-file operations
	
	
					
	//========================================================================== two-file operations
	if (mode == "add" || mode == "sub" || mode == "mult" || mode == "div"
			|| mode == "mask_apply" ){ 
		if (num_files > 1){
			int fail = 0;
			fail += io->readFromFile( in_fn.at(0), first );
			fail += io->readFromFile( in_fn.at(1), second );
			if (fail){
				cerr << "Error reading file" << endl;
				return 3;
			}
		}else{
			cerr << "Not enough input files for this mode" << endl;
			return 3;
		}
		
		if (mode == "add"){
			cout << "ADDITION: " << in_fn.at(0) << " + " << in_fn.at(1) << "" << endl;
			first->addArrayElementwise( second );
		}
		
		if (mode == "sub"){
			cout << "SUBTRACTION: " << in_fn.at(0) << " - " << in_fn.at(1) << "" << endl;
			first->subtractArrayElementwise( second );
		}
		
		if (mode == "mult"){
			cout << "MULTIPLICATION: " << in_fn.at(0) << " * " << in_fn.at(1) << "" << endl;
			first->multiplyByArrayElementwise( second );
		}
		
		if (mode == "div"){
			cout << "DIVISION: " << in_fn.at(0) << " / " << in_fn.at(1) << "" << endl;
			first->divideByArrayElementwise( second );
		}
		
		if (mode == "mask_apply"){
			cout << "APPLICATION OF A MASK" << endl;
			first->applyMask( second );
		}
		
		//write output
		io->writeToFile( out_fn, first );
	}//end of two-file operations

	
	//========================================================================== any-number operations
	if (mode == "stats" || mode == "avg" || mode == "IO"){
	
		// 1) do some general setup before entering the loop
		array2D *sum = 0;
		if (mode == "stats"){
			io->setVerbose(0);
			//write header
			std::ostringstream info_global;	
			info_global << "no "  
				<< std::setw(10) << "min" 
				<< std::setw(10) << "max" 
				<< std::setw(10) << "mean" 
				<< std::setw(10) << "stdev" 
				<< std::setw(10) << "pix1" 
				<< std::setw(10) << "pix2"
				<< std::setw(10) << "pix3"
				<< " filename "
				;
			cout << info_global.str() << endl;
			fout << info_global.str() << endl;
		}
		
		// 2) main loop: go through all files
		for (int i = 0; i < num_files; i++){
			std::ostringstream info;	
			array2D *data = 0;
			int fail = io->readFromFile( in_fn.at(i), data );
			if (fail){
				cerr << "Error reading file" << endl;
				return 3;
			}
			
			if (mode == "avg"){
				if (i == 0){
					sum = new array2D(*data);
				}else{
					sum->addArrayElementwise(data);
				}
			}
			
			if (mode == "stats"){
				double min = data->calcMin();
				double max = data->calcMax();
				double mean = 0;
				double stdev = data->calcStDev(mean);
				double pixel1 = data->get_atIndex(1000);
				double pixel2 = data->get_atIndex(2000);
				double pixel3 = data->get_atIndex(3000);
				
				info << i
					<< std::setw(10) << min
					<< std::setw(10) << max
					<< std::setw(10) << mean
					<< std::setw(10) << stdev
					<< std::setw(10) << pixel1					
					<< std::setw(10) << pixel2
					<< std::setw(10) << pixel3
					<< " " << in_fn.at(i)
					;
					
				cout << info.str() << endl;
				fout << info.str() << endl;
			}
			
			if (mode == "IO"){
				std::ostringstream osst;
				osst << in_fn.at(i) << "_" << out_fn;
				io->writeToFile( osst.str(), data );
			}
			
			delete data;
		}
		
		// 3) do some post-processing
		if (mode == "avg"){
			sum->divideByValue(num_files);
			io->writeToFile( out_fn, sum );
			delete sum;
		}
		
	}//end of any-number file operations

	//clean up
	fout.close();
	delete first;
	delete second;
	delete io;

	return 0;
}

