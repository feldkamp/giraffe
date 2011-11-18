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
    //cout << "Hello, World!\n";

	int alg = 1;
	string list_fn = "";
	string file_fn = "";
	string pixx_fn = "";
	string pixy_fn = "";
	string outdir = "";
	string back_fn = "";
	string mask_fn = "";
	bool singleout = false;
	double weight = 0.;

	po::options_description desc("Allowed options");
	desc.add_options()		
		("help,h",																"produce help message")		
		("file,f", 			po::value<string>(&file_fn), 						"single primary input file")
		("list,l", 			po::value<string>(&list_fn), 						"file list with input files")
		("pixelX", 			po::value<string>(&pixx_fn), 						"input file for pixel values x-coordinate")
		("pixelY", 			po::value<string>(&pixy_fn), 						"input file for pixel values y-coordinate")
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
		po::notify(vm);
	}catch( const boost::program_options::error &e ){
		cerr << "Error in command line arguments: " << e.what() << "." << endl;
		cerr << "Use option '--help' or '-h' to see a list of valid options." << endl;
		exit(1);
	}	
	//cout << ">parser done<" << endl;
	
	if (vm.count("help")) {
		cout << desc << endl;
		exit(0);
	}
	if (vm.count("f")){
		cout << "--> using input file " << file_fn << endl;
	}
	if (vm.count("pixelX")){
		cout << "--> using x-values file " << pixx_fn << endl;
	}
	if (vm.count("pixelY")){
		cout << "--> using y-values file " << pixy_fn << endl;
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
	
	//get primary input file(s) and store them in 'files'
	std::vector<string> files;
		
	if (list_fn == "" && file_fn == ""){
		cerr << "Error. No input file(s) given. Use -f or -l to specify a (f)ile or a (l)ist of files" << endl;
		exit(2);
	}else if (list_fn != ""){
		cout << "--> using file list in '" << list_fn << "'" << endl;
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
			cerr << "Error. Could not open file list '" << list_fn << "', aborting..." << endl;
			exit( 1 );
		}
		fin.close();
	}else{
		cout << "--> using file '" << file_fn << "'" << endl;
		files.push_back(file_fn);
	}
		
	//create analyzer object with the given options
	Analyzer *ana = new Analyzer();
	ana->setAlg( alg );	
	ana->setOutputDirectory( outdir );
	ana->setFlagSingleCorrelationOutput( true );
	ana->setBackgroundWeight( weight );
	
	//define set of variables to pass to the processFiles function, should probabaly go into an ini file or so at some point
	int num_phi = 512;
	int num_q = 400;
	double start_q = 0;			//units of these should match those in pixx, pixy!
	double stop_q = 800;
	int LUTx = 1000;
	int LUTy = 1000;

	int fail = ana->processFiles( files, num_phi, num_q, start_q, stop_q, LUTx, LUTy, 
						pixx_fn, pixy_fn, mask_fn, back_fn);

    delete ana;
    return fail;
}



