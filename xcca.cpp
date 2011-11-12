#include "analyze.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
using std::flush;

#include <string>
using std::string;

#include <sstream>		//strings
#include <fstream>		//files
#include <vector>

#include <boost/program_options.hpp>		//parser for options (commandline and ini-file)
namespace po = boost::program_options;

#include "crosscorrelator.h"
#include "arrayclasses.h"
#include "arraydataIO.h"


int main (int argc, char * const argv[]) {
    cout << "Hello, World!\n";

	int alg = 0;
	string list_fn = "";
	string outdir = "";
	string back_fn = "";
	string mask_fn = "";
	bool singleout = false;
	double weight = 0.;

	po::options_description desc("Allowed options");
	desc.add_options()		
		("help,h",																"produce help message")		
		("list,l", 			po::value<string>(&list_fn), 						"file list with input files")
		("alg,a", 			po::value<int>(&alg),								"set algorithm")
		("outdir,o", 		po::value<string>(&outdir),							"output directory")
		("back,b", 			po::value<string>(&back_fn), 						"background file")
		("mask,m", 			po::value<string>(&mask_fn), 						"mask file")
		("single,s", 		po::bool_switch(&singleout), 						"single image output?")
		("weight,B", 		po::value<double>(&weight),					 		"background weighting factor")
	;
	
	po::variables_map vm;
	try{
		po::parsed_options parsed = po::parse_command_line(argc, argv, desc);
		//po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
		po::store(parsed, vm);
	}catch(...){
		cerr << "Exception!" << endl;
	}
	po::notify(vm);
	
	cout << ">parser done<" << endl;
	
	if (vm.count("help")) {
		cout << desc << endl;
	}
	if (vm.count("a")) {
		//cout << "algorithm " << vm["alg"].as<int>() << "." << endl;
		cout << "--> using algorithm " << alg << endl;
	}
	if (vm.count("s")){
		cout << "--> using single image output " << endl;
	}
	if (vm.count("o")){
		cout << "--> using output directory '" << outdir << "'" << endl;
	}
	if (vm.count("b")){
		cout << "--> using background file '" << back_fn << "'" << endl;
	}
	if (vm.count("m")){
		cout << "--> using mask file '" << mask_fn << "'" << endl;
	}
	if (vm.count("B")){
		cout << "--> using background weighting factor " << weight << endl;	
	}
	


	if (list_fn == ""){
		cerr << "no file list given. use -l to specify a file list" << endl;
		exit(2);
	}else{
		cout << "--> using file list in '" << list_fn << "'" << endl;
	}
	
	//read from file list
	std::vector<string> files;
	std::ifstream fin;
	fin.open(list_fn.c_str());					
	if( !fin.fail() ){
		string line;
		while (fin.good()) {
			getline(fin, line);
			if (line != ""){
				cout << line << endl;
				files.push_back(line);
			}
		}
	}else{
		cerr << "could not open file list '" << list_fn << "', aborting..." << endl;
		exit( 1 );
	}
	fin.close();


	arraydataIO *io = new arraydataIO;

		
	//create analyzer object with the given options
	Analyzer *ana = new Analyzer();
	ana->setAlg( alg );	
	ana->setOutputDirectory( outdir );
	ana->setFlagSingleCorrelationOutput( true );
	ana->setBackgroundWeight( weight );
	
	if (back_fn != ""){
		array2D *back = new array2D;
		io->readFromFile( back_fn, back);
		ana->setBackground( back );
		ana->setFlagSubtractBackground( true );	
		delete back;
	}
	
	cout << "." << flush;
	
	if (mask_fn != ""){
		array2D *mask = new array2D;
		io->readFromFile( mask_fn, mask);
		ana->setMask( mask );			
		delete mask;
	}
	
	cout << "." << flush;
	
	//define set of variables to pass to the processFiles function, should probabaly go into an ini file or so at some point
	double shiftX = -7;
	double shiftY = 1;
	int num_phi = 512;
	int num_q = 200;
	double start_q = 0;
	double stop_q = 800;
	int LUTx = 487;				//pilatus detector x and y values
	int LUTy = 619;


	int fail = ana->processFiles( files, shiftX, shiftY, num_phi, num_q, start_q, stop_q, LUTx, LUTy);

	cout << "Goodbye, " << flush;

    delete ana;
	delete io;
	
    cout << "World" << endl;
    return fail;
}



