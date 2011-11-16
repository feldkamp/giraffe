//
//  test.cpp
//  xcca_commandline
//
//  Created by Feldkamp on 7/21/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "analyze.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

#include <sstream>
#include <string>
using std::string;

#include <vector>
using std::vector;

#include "crosscorrelator.h"

Analyzer::Analyzer()
	: io()
	, p_back_weight(1.)
	, p_out_dir("")
	, p_alg(1)
	, p_flag_subtract_background(false)
	, p_flag_single_correlation_output(false)
	, p_flag_use_mask(false)
{
	io = new arraydataIO;
}

Analyzer::~Analyzer(){
	delete io;
}

// =============================================================
// assuming all files are similar
// i.e. same dimensions, same beam center
//
//	int number_phi = 256;
//	int number_q = 200;
//	int cenX = 260;			//in detector pixels
//	int cenY = 308;	
// =============================================================
int Analyzer::processFiles( vector<string> files, int num_phi, int num_q){
	double start_q = 0;
	double stop_q = num_q;
	int LUTx = 200;
	int LUTy = 200;
	return processFiles(files, num_phi, num_q, start_q, stop_q, LUTx, LUTy, 
							"", "", "", "");
}

int Analyzer::processFiles( vector<string> files, int num_phi, int num_q, double start_q, double stop_q, int LUTx, int LUTy, 
							string qx_fn, string qy_fn, string mask_fn, string back_fn ){
	
	array2D *mask = 0;
	array2D *back = 0;
	array2D *qx = 0;
	array2D *qy = 0;
	
	int nQ = 0;
	int nPhi = 0;
	int nLag = 0;

	double shiftX = 0.;
	double shiftY = 0.;
	
	//pilatus specifics for P06 July 2011
	//	shiftX = -7;
	//	shiftY = 1;


	//prepare array2D's to hold overall averages
	array2D *detavg = new array2D;
	io->readFromFile( files.at(0), detavg );
	int imgX = detavg->dim2();
	int imgY = detavg->dim1();
	
	array2D *backavg = new array2D( imgY, imgX );
	
			
	if (mask_fn != ""){
		io->readFromFile( mask_fn, mask);
		setFlagUseMask( true );
	}
	
	if (back_fn != ""){
		io->readFromFile( back_fn, back);
		setFlagSubtractBackground( true );	
	}
	
	if (qx_fn != ""){
		io->readFromFile( qx_fn, qx);
	}else{
		//generate a qx-distribution
		qx = new array2D( imgY, imgX );
		qx->gradientAlongDim1(-imgX/2+shiftX, +imgX/2+shiftX);
	}

	if (qy_fn != ""){
		io->readFromFile( qy_fn, qy);
	}else{
		//generate a qy-distribution
		qy = new array2D( imgY, imgX );
		qy->gradientAlongDim2(-imgY/2+shiftY, +imgY/2+shiftY);
	}
	
	
	string ext = ".h5";			//standard extension for output (alternatives: .edf, .tif, .txt)
	string outDir = this->outputDirectory();
	cout << "Output will be written to " << (outDir=="" ? "current directory" : outDir) 
		<< " with extension '" << ext << "'" << endl;

	if (files.size() == 0){
		cerr << "no files found. aborting." << endl;
		return 1;
	}


	//prepare lookup table once, so it doesn't have to be done every time
	CrossCorrelator *dummy_cc= new CrossCorrelator(detavg, qx, qy, num_phi, num_q);
	dummy_cc->createLookupTable(LUTy, LUTx);
	array2D *LUT = new array2D( *(dummy_cc->lookupTable()) );
	//io->writeToFile( outDir+"LUT"+ext, LUT );
	nQ = dummy_cc->nQ();
	nPhi = dummy_cc->nPhi();
	nLag = dummy_cc->nLag();
	delete dummy_cc;
	array2D *polaravg = new array2D( nQ, nPhi );
	array2D *corravg = new array2D( nQ, nLag );
	
	if ( flagSubtractBackground() ){
		//prepare background file
		back->multiplyByValue( backgroundWeight() );
	}
	
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
		
		
		CrossCorrelator *cc = new CrossCorrelator(image, qx, qy, num_phi, num_q);
		if ( flagUseMask() ){
			cc->setMask( mask );
		}
		cc->setOutputdir( outDir );
		cc->setDebug(0);
		
		if (flagSubtractBackground()){
			image->subtractArrayElementwise( back );
			backavg->addArrayElementwise( back );
		}
		
		//run the calculation...
		bool calcSAXS = true;
		cc->run(start_q, stop_q, alg(), calcSAXS);

		if ( flagSingleCorrelationOutput() ){
			io->writeToFile( outDir+"corr"+single_desc+ext, cc->autoCorr() );
			io->writeToFile( outDir+"polar"+single_desc+ext, cc->polar() );
			io->writeToFile( outDir+"pixel_count"+single_desc+ext, cc->pixelCount() );
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
	
	io->writeToFile( outDir+"det_avg"+ext, detavg);			// average background-subtracted detector image
	io->writeToFile( outDir+"polar_avg"+ext, polaravg);		// average image in polar coordinates
	io->writeToFile( outDir+"corr_avg"+ext, corravg);			// average autocorrelation
	if ( flagSubtractBackground() ){
		io->writeToFile( outDir+"det_background_avg"+ext, backavg);					// just the background
	}

	delete detavg;
	delete backavg;
	delete polaravg;
	delete corravg;		
	delete LUT;
	
	return 0;
}


void Analyzer::setBackgroundWeight( double weight ){
	p_back_weight = weight;
}

double Analyzer::backgroundWeight(){
	return p_back_weight;
}


void Analyzer::setOutputDirectory( string outdir ){
	if (outdir != ""){
		//if last character is not a slash, append one
		const char lastchar = outdir.at( outdir.size()-1 );
		if( lastchar != '/' ){		
			outdir += '/';
		}
	}
	p_out_dir = outdir;
}

string Analyzer::outputDirectory(){
	return p_out_dir;
}

void Analyzer::setAlg( int alg ){
	p_alg = alg;
}

int Analyzer::alg(){
	return p_alg;
}

void Analyzer::setFlagSubtractBackground( bool flag ){
	p_flag_subtract_background = flag;
}

bool Analyzer::flagSubtractBackground(){
	return p_flag_subtract_background;
}

void Analyzer::setFlagSingleCorrelationOutput( bool flag ){
	p_flag_single_correlation_output = flag;
}

bool Analyzer::flagSingleCorrelationOutput(){
	return p_flag_single_correlation_output;
}

void Analyzer::setFlagUseMask( bool flag ){
	p_flag_use_mask = flag;
}

bool Analyzer::flagUseMask(){
	return p_flag_use_mask;
}

