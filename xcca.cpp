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
	bool single_out = false;
	double back_weight = 0.;
	
	int nPhi = 512;
	int nLag = nPhi/2. + 1;
	int nQ = 400;
	int LUTx = 1000;
	int LUTy = 1000;
	
	double start_q = 0;			//units of these should match those in pixx, pixy!
	double stop_q = nQ;

	double shiftX = 0.;
	double shiftY = 0.;
	
	//pilatus specifics for P06 July 2011
	//	shiftX = -7;
	//	shiftY = 1;

	po::options_description desc("Allowed options");
	desc.add_options()		
		("help,h",																"produce help message")		
		("file,f", 			po::value<string>(&file_fn), 						"single primary input file")
		("list,l", 			po::value<string>(&list_fn), 						"file list with input files")
		("pixelX", 			po::value<string>(&pixx_fn), 						"input file for pixel values x-coordinate")
		("pixelY", 			po::value<string>(&pixy_fn), 						"input file for pixel values y-coordinate")
		("numPhi", 			po::value<int>(&nPhi),								"number of angles phi")
		("numQ", 			po::value<int>(&nQ),								"number of q-values")
		("alg,a", 			po::value<int>(&alg),								"set algorithm")
		("outdir,o", 		po::value<string>(&outdir),							"output directory")
		("back,b", 			po::value<string>(&back_fn), 						"background file")
		("mask,m", 			po::value<string>(&mask_fn), 						"mask file")
		("single,s", 		po::bool_switch(&single_out), 						"single image output?")
		("weight,B", 		po::value<double>(&back_weight),					"background weighting factor")
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
		cout << "--> using background weighting factor " << back_weight << endl;	
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

	//-------------------------------------------------interface done, get to work

	arraydataIO *io = new arraydataIO;

	//prepare array2D's to hold overall averages
	array2D *detavg = new array2D;
	io->readFromFile( files.at(0), detavg );
	int imgX = detavg->dim2();
	int imgY = detavg->dim1();
	
	array2D *backavg = new array2D( imgY, imgX );
	array2D *mask = 0;
	array2D *back = 0;
	array2D *qx = 0;
	array2D *qy = 0;
				
	if (mask_fn != ""){
		io->readFromFile( mask_fn, mask);
	}
	
	if (back_fn != ""){
		io->readFromFile( back_fn, back);
		if ( back_weight != 1 ){
			back->multiplyByValue( back_weight );
		}
	}
	
	if (pixx_fn != ""){
		io->readFromFile( pixx_fn, qx);
	}else{
		//generate a qx-distribution
		qx = new array2D( imgY, imgX );
		qx->gradientAlongDim1(-imgX/2+shiftX, +imgX/2+shiftX);
	}

	if (pixy_fn != ""){
		io->readFromFile( pixy_fn, qy);
	}else{
		//generate a qy-distribution
		qy = new array2D( imgY, imgX );
		qy->gradientAlongDim2(-imgY/2+shiftY, +imgY/2+shiftY);
	}
	
	
	string ext = ".h5";			//standard extension for output (alternatives: .edf, .tif, .txt)
	cout << "Output will be written to " << (outdir=="" ? "current directory" : outdir) 
		<< " with extension '" << ext << "'" << endl;

	if (files.size() == 0){
		cerr << "no files found. aborting." << endl;
		return 1;
	}


	//prepare lookup table once, so it doesn't have to be done every time
	CrossCorrelator *dummy_cc= new CrossCorrelator(detavg, qx, qy, nPhi, nQ);
	dummy_cc->createLookupTable(LUTy, LUTx);
	array2D *LUT = new array2D( *(dummy_cc->lookupTable()) );
	//io->writeToFile( outDir+"LUT"+ext, LUT );
	nLag = dummy_cc->nLag();
	delete dummy_cc;
	array2D *polaravg = new array2D( nQ, nPhi );
	array2D *corravg = new array2D( nQ, nLag );
	
	//process all files
	unsigned int num_files = (unsigned int)files.size();
	for (int k = 0; k < num_files; k++){
		std::ostringstream osst_num;
		osst_num << k;
		string single_desc = osst_num.str();
		string fn = files.at(k);
		cout << "#" << k << ": " << std::flush;// << fn << endl;
		
		array2D *image = 0;
		if (k != 0){
			io->readFromFile( fn, image );
		}else{
			//this has been read previously in detavg, no need to read again...
			delete image;
			image = new array2D(*detavg);
			detavg->zeros();
		}
		
		
		CrossCorrelator *cc = new CrossCorrelator(image, qx, qy, nPhi, nQ);
		if ( mask ){
			cc->setMask( mask );
		}
		cc->setOutputdir( outdir );
		cc->setDebug(0);
		
		if (back){
			image->subtractArrayElementwise( back );
			backavg->addArrayElementwise( back );
		}
		
		//run the calculation...
		bool calcSAXS = true;
		cc->run(start_q, stop_q, alg, calcSAXS);

		if ( single_out ){
			io->writeToFile( outdir+"corr"+single_desc+ext, cc->autoCorr() );
			io->writeToFile( outdir+"polar"+single_desc+ext, cc->polar() );
			io->writeToFile( outdir+"pixel_count"+single_desc+ext, cc->pixelCount() );
		}
		
		//sum up results
		polaravg->addArrayElementwise( cc->polar() );
		corravg->addArrayElementwise( cc->autoCorr() );
		detavg->addArrayElementwise( image );
		
		delete image;
		delete cc;
	}
	
	detavg->divideByValue( num_files );			//normalize
	backavg->divideByValue( num_files );		//normalize
	polaravg->divideByValue( num_files );		//normalize	
	corravg->divideByValue( num_files );		//normalize
	
	io->writeToFile( outdir+"det_avg"+ext, detavg);			// average background-subtracted detector image
	io->writeToFile( outdir+"polar_avg"+ext, polaravg);		// average image in polar coordinates
	io->writeToFile( outdir+"corr_avg"+ext, corravg);			// average autocorrelation
	if ( back ){
		io->writeToFile( outdir+"det_background_avg"+ext, backavg);		// just the background
	}

	delete detavg;
	delete backavg;
	delete polaravg;
	delete corravg;		
	delete LUT;	
	delete io;
	
    return 0;
}



