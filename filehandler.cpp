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
#include "util.h"




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
	std::ostringstream info;					// contents of info string are written to logfile
	
	string mode = "";	
	vector<string> in_fn;
	string list_fn = "";
	string out_fn = "";
	double p1 = 0.;
	double p2 = 0.;
	int verbose = 0;

	po::positional_options_description p;
	p.add("mode", -1);

	po::options_description desc("Usage: filehandler <mode> [options] \nAllowed options:");
	desc.add_options()		
		("mode,m", 		po::value<string>(&mode),								"mode of operation,"
									"\nchoose a keyword:"
									"\n   add                                       (2 files)"
									"\n   sub                                       (2 files)"
									"\n   mult                                      (2 files)"
									"\n   div                                       (2 files)"
									"\n   thresholding --param_one --param_two      (1 file)"
									"\n       param_one = lower threshold, param_two = upper threshold"
									"\n   mask_invert                               (1 file)"
									"\n   mask_apply                                (2 files)"
									"\n   mask_combine                              (x files)"
									"\n   stats                                     (x files)"
									"\n   avg                                       (x files)"
									"\n   IO                                        (x files)"
									"\n      output file type is determined by extension"
									"\n   cspad_raw_img                             (x files)"
									)
		("param_one", 	po::value<double>(&p1), 								"generic parameter 1 (use depends on chosen mode)")
		("param_two", 	po::value<double>(&p2), 								"generic parameter 2 (use depends on chosen mode)")
		("output,o",	po::value<string>(&out_fn), 							"output filename")	
		("input,i", 	po::value< vector<string> >(&in_fn)->multitoken(),		"input filename")
		("list,l", 		po::value<string>(&list_fn),							"input filelist")
		("verbose,v", 	po::value<int>(&verbose), 								"set verbosity level")
		("help,h",																"display help message")	
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

	if (vm.count("verbose")){
		cout << "--> using verbosity " << verbose << endl;
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
	if (mode == "thresholding" || mode == "mask_invert"){ 
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

		if (mode == "thresholding"){
			cout << "THRESHOLDING" << endl;
			double lower = p1;
			double upper = p2;
			int keep, use_high, use_low = 0;
			array2D *rejected = new array2D( first->dim1(), first->dim2() );
			for (int i=0; i<first->size(); i++){
				if (first->get_atIndex(i) > lower){
					if (first->get_atIndex(i) <= upper){
						//keep original value
						rejected->set_atIndex(i, 1);
						keep++;
					} else {
						//use upper value
						first->set_atIndex(i, upper);
						rejected->set_atIndex(i, 0);
						use_high++;
					}
				} else {
					//use lower value
					first->set_atIndex(i, lower);
					rejected->set_atIndex(i, 0);
					use_low++;
				}
			}
			string rejected_fn = "rejected_px_after_thresholding.h5";
			io->writeToHDF5( rejected_fn, rejected );
			delete rejected;
			double total = (double)first->size();
			cout << total << " values: " << keep << " kept, " 
				<< use_high << " (" << use_high/total*100. << "%) replaced by high limit " << upper << ", "
				<< use_low << " (" << use_low/total*100. << "%) replaced by low limit " << lower << "." << endl;
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
	if (mode == "stats" || mode == "avg" || mode == "IO" || mode == "mask_combine" || mode == "cspad_rawimg"){
	
		// 1) do some general setup before entering the loop
		array2D *sum = 0;
		int dim1 = 0;
		int dim2 = 0;
		if (mode == "stats"){
			std::ostringstream header;	
			io->setVerbose(0);

			header << "no "  
				<< std::setw(10) << "min" 
				<< std::setw(10) << "max" 
				<< std::setw(10) << "mean" 
				<< std::setw(10) << "stdev" 
				<< std::setw(10) << "pix1" 
				<< std::setw(10) << "pix2"
				<< std::setw(10) << "pix3"
				<< " filename "
				;
				
			cout << header.str() << endl;
			info << header.str() << endl;
		}
		
		// 2) main loop: go through all files
		for (int i = 0; i < num_files; i++){
			std::ostringstream osst;
			osst << in_fn.at(i) << "_" << out_fn;
			string current_out_fn = osst.str();
			if (num_files == 1){
				current_out_fn = out_fn;		// if there's only one file, don't change the output file name
			}
			
			std::ostringstream stat;	
			array2D *data = 0;
			int fail = io->readFromFile( in_fn.at(i), data );
			if (fail){
				cerr << "Error reading file" << endl;
				return 3;
			}
			if (i == 0){
				dim1 = data->dim1();
				dim2 = data->dim2();
			}else if (dim1 != data->dim1() || dim2 != data->dim2()) {
				cerr << "WARNING in filehandler: Dimensions of all input files don't match! " 
					<< "First file: " << dim1 << "," << dim2 
					<< ", this file: " << data->dim1() << "," << data->dim2() << endl;
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
				
				stat << i
					<< std::setw(10) << min
					<< std::setw(10) << max
					<< std::setw(10) << mean
					<< std::setw(10) << stdev
					<< std::setw(10) << pixel1					
					<< std::setw(10) << pixel2
					<< std::setw(10) << pixel3
					<< " " << in_fn.at(i)
					;
					
				cout << stat.str() << endl;
				info << stat.str() << endl;
			}
			
			if (mode == "IO"){
				io->writeToFile( current_out_fn, data );
			}
			
			if (mode == "mask_combine"){
				//multiple execution of applyMask
				if (i == 0){
					sum = new array2D(*data);
				}else{
					sum->applyMask( data );				
				}
			}
			
			if (mode == "cspad_rawimg"){
				//transpose() causes to essentially read rows in the file first
				//this is the way it needs to be, for example, with the files generated by the pedestal module
				data->transpose();		
				array2D *img = 0;
				ns_cspad_util::createRawImageCSPAD( data, img );
				io->writeToFile( current_out_fn, img );
				delete img;
			}
			
			delete data;
		}
		
		// 3) do some post-processing
		if (mode == "avg"){
			sum->divideByValue(num_files);
			io->writeToFile( out_fn, sum );
			delete sum;
		}
		
		if (mode == "mask_combine"){
			io->writeToFile( out_fn, sum );
		}
		
	}//end of any-number file operations

	if (info.str() != ""){
		std::ofstream fout("filehandler.log");
		fout << info.str() << endl;	
		fout.close();
	}
	
	//clean up
	delete first;
	delete second;
	delete io;

	return 0;
}

