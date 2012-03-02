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

#include "arrayclasses.h"
using namespace ns_arraydata;
#include "crosscorrelator.h"
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
	int verbose = 0;
	
	int nPhi = 512;
	int nLag = (int)(nPhi/2. + 1);
	int nQ = 100;
	
	//units of start_q, stop_q should match those in pixx, pixy!
	double start_q = 0;			
	double stop_q = 0;

	double shiftX = 0.;
	double shiftY = 0.;
	
	//pilatus specifics for P06 July 2011
	//	shiftX = -7;
	//	shiftY = 1;

	po::options_description desc("Allowed options");
	desc.add_options()		
		("help,h",																"produce help message")		
		("input,i", 		po::value<string>(&file_fn), 						"single primary input file")
		("list,l", 			po::value<string>(&list_fn), 						"file list with input files")
		("pixelX", 			po::value<string>(&pixx_fn), 						"input file for pixel values x-coordinate")
		("pixelY", 			po::value<string>(&pixy_fn), 						"input file for pixel values y-coordinate")
		("shiftX", 			po::value<double>(&shiftX), 						"number of pixels in x to shift image center by (when not given pixelX)")
		("shiftY", 			po::value<double>(&shiftY), 						"number of pixels in y to shift image center by (when not given pixelY)")
		("numPhi", 			po::value<int>(&nPhi),								"number of angles phi")
		("numQ", 			po::value<int>(&nQ),								"number of q-values")
		("startQ", 			po::value<double>(&start_q),							"start value for q")
		("stopQ", 			po::value<double>(&stop_q),							"stop value for q")
		("alg,a", 			po::value<int>(&alg),								"set algorithm")
		("outdir,o", 		po::value<string>(&outdir),							"output directory")
		("back,b", 			po::value<string>(&back_fn), 						"background file")
		("mask,m", 			po::value<string>(&mask_fn), 						"mask file")
		("single,s", 		po::bool_switch(&single_out), 						"single image output?")
		("weight,B", 		po::value<double>(&back_weight),					"background weighting factor")
		("verbose,v", 		po::value<int>(&verbose), 							"set verbosity level")
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
	if (vm.count("input")){
		cout << "--> using input file " << file_fn << endl;
	}
	if (vm.count("pixelX")){
		cout << "--> using x-values file " << pixx_fn << endl;
	}
	if (vm.count("pixelY")){
		cout << "--> using y-values file " << pixy_fn << endl;
	}
	if (vm.count("shiftX")){
		cout << "--> shifting image center in x by " << shiftX << endl;
	}
	if (vm.count("shiftY")){
		cout << "--> shifting image center in y by " << shiftY << endl;
	}
	if (vm.count("numPhi")){
		cout << "--> using numPhi value " << nPhi << endl;
	}
	if (vm.count("numQ")){
		cout << "--> using numQ value " << nQ << endl;
	}
	if (vm.count("startQ")){
		cout << "--> using startQ value " << start_q << endl;
	}
	if (vm.count("stopQ")){
		cout << "--> using stopQ value " << stop_q << endl;
	}
	if (vm.count("a")) {
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
	if (vm.count("verbose")){
		cout << "--> using verbosity " << verbose << endl;
	}	
	
	//get primary input file(s) and store them in 'files'
	std::vector<string> files;
		
	if (list_fn == "" && file_fn == ""){
		cerr << "Error. No input file(s) given. Use -f or -l to specify a (f)ile or a (l)ist of files" << endl;
		cout << desc << endl;
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
	arraydataIO *io = new arraydataIO(verbose);
	int iofail = 0;
	
	//prepare array2D<double>'s to hold overall averages
	array2D<double> *detavg = new array2D<double>;
	iofail += io->readFromFile( files.at(0), detavg );
	int imgX = detavg->dim2();
	int imgY = detavg->dim1();
	
	array2D<double> *backavg = new array2D<double>( imgY, imgX );
	array2D<double> *mask = 0;
	array2D<double> *back = 0;
	array2D<double> *qx = 0;
	array2D<double> *qy = 0;
	

				
	if (mask_fn != ""){
		iofail += io->readFromFile( mask_fn, mask);
	}
	
	if (back_fn != ""){
		iofail += io->readFromFile( back_fn, back);
		if ( back_weight != 1 ){
			back->multiplyByValue( back_weight );
		}
	}
	
	if (pixx_fn != ""){
		iofail += io->readFromFile( pixx_fn, qx);
	}else{
		cout << "using a generic qx-distribution (with " << shiftX << " shift)" << endl;
		qx = new array2D<double>( imgY, imgX );
		qx->gradientAlongDim1(-imgX/2+shiftX, +imgX/2+shiftX);
	}

	if (pixy_fn != ""){
		iofail += io->readFromFile( pixy_fn, qy);
	}else{
		cout << "using a generic qy-distribution (with " << shiftY << " shift)" << endl;
		qy = new array2D<double>( imgY, imgX );
		qy->gradientAlongDim2(-imgY/2+shiftY, +imgY/2+shiftY);
	}
	
	if (iofail){
		cerr << "I/O failure. Aborting." << endl;
		return 1;
	}
	
	//if stop_q is still zero, give it a reasonable default value
	if ( stop_q == 0){
		stop_q = qx->calcMax();
	}
	
	string ext = ".h5";			//standard extension for output (alternatives: .edf, .tif, .txt)
	cout << "Output will be written to " << (outdir=="" ? "current directory" : outdir) 
		<< " with extension '" << ext << "'" << endl;

	if (files.size() == 0){
		cerr << "no input files found. aborting." << endl;
		exit(2);
	}


	//prepare lookup table once, so it doesn't have to be done every time
	CrossCorrelator *cc= new CrossCorrelator(detavg, qx, qy, nPhi, nQ);
	
	if ( mask ){
		cc->setMask( mask );
	}
	cc->setOutputdir( outdir );
	cc->setDebug(verbose);
		
	if (alg == 2 || alg == 4) {
		int LUTx = 1000;
		int LUTy = 1000;
		cc->createLookupTable(LUTy, LUTx);
	}

	nLag = cc->nLag();
	
	array2D<double> *polaravg = new array2D<double>( nQ, nPhi );
	array2D<double> *corravg = new array2D<double>( nQ, nLag );
	
	//process all files
	unsigned int num_files = (unsigned int)files.size();
	for (int k = 0; k < num_files; k++){
		std::ostringstream osst_num;
		osst_num << k;
		string single_desc = "evt_"+osst_num.str();
		string fn = files.at(k);
		cout << "#" << k << ": " << std::flush;// << fn << endl;
		
		array2D<double> *image = 0;
		if (k != 0){
			iofail = io->readFromFile( fn, image );
		}else{
			//this has been read previously in detavg, no need to read again...
			delete image;
			image = new array2D<double>(*detavg);
			detavg->zeros();
		}
		
		if(iofail){
			cerr << "I/O failure. Aborting." << endl;
			return 1;		
		}
		
		
		if (back){
			image->subtractArrayElementwise( back );
			backavg->addArrayElementwise( back );
		}
		
		//run the calculation...
		try{
			cc->setData(image);
			cc->run(start_q, stop_q, alg);
		}
		catch(...){
			cerr << "\nException caught while running the cross correlation. Aborting..." << endl;
			exit(1);
		}

		if ( single_out || num_files == 1 ){
			io->writeToFile( outdir+single_desc+"xaca"+ext, cc->autoCorr() );
			io->writeToFile( outdir+single_desc+"polar"+ext, cc->polar() );
			io->writeToFile( outdir+single_desc+"fluct"+ext, cc->fluctuations() );
			
			//pixelCount is of type array2D<unsigned int> and would need to be converted to be written
			//io->writeToFile( outdir+single_desc+"pixel_count"+ext, cc->pixelCount() );
			
			io->writeToFile( outdir+single_desc+"q"+ext, cc->qAvg() );
			io->writeToFile( outdir+single_desc+"i"+ext, cc->iAvg() );
			
			io->writeToFile( outdir+single_desc+"det"+ext, cc->data() );
		}
		
		//sum up results
		if (num_files > 1){
			polaravg->addArrayElementwise( cc->polar() );
			corravg->addArrayElementwise( cc->autoCorr() );
			detavg->addArrayElementwise( image );
		}
			
		delete image;
	}
	
	
	//compute and output average if there was more than one file
	if (num_files > 1){
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
	}
	
	delete detavg;
	delete backavg;
	delete polaravg;
	delete corravg;		
	delete io;
	delete cc;
	
    return 0;
}



