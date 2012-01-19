//
//  util.cpp
//  ana
//
//  Created by Feldkamp on 9/28/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "util.h"
using ns_cspad_util::nRowsPerASIC;
using ns_cspad_util::nColsPerASIC;
using ns_cspad_util::nRowsPer2x1;
using ns_cspad_util::nColsPer2x1;
using ns_cspad_util::nMax2x1sPerQuad;
using ns_cspad_util::nMaxQuads;
using ns_cspad_util::nPxPer2x1;
using ns_cspad_util::nMaxPxPerQuad;
using ns_cspad_util::nMaxTotalPx;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <string>
using std::string;

#include <sstream>

#include <cmath>


namespace ns_cspad_util{

	//------------------------------------------------------------- create1DFromRawImageCSPAD
	// inverse function of "createRawImageCSPAD"
	// make a 1D array from a written 2D "raw" image
	//
	//-------------------------------------------------------------
	int create1DFromRawImageCSPAD( const array2D<double> *input, array1D<double> *&output ){
		if (!input){
			cerr << "Error in create1DFromRawImageCSPAD. No input. Aborting..." << endl;
			return 1;
		}
		
		delete output;
		output = new array1D<double>( nMaxTotalPx );

		//transpose to be conform with psana's convention
		//make a copy to not screw up the original array
		//obviously not the smart way to do it... needs a fix at a certain point...... should be simple by rearranging logic in 4-for-loop
		array2D<double> *copy = new array2D<double>( *input );
		copy->transpose();

		//sort into 1D data
		for (int q = 0; q < nMaxQuads; q++){
			int superrow = q*nRowsPer2x1;
			for (int s = 0; s < nMax2x1sPerQuad; s++){
				int supercol = s*nColsPer2x1;
				for (int c = 0; c < nColsPer2x1; c++){
					for (int r = 0; r < nRowsPer2x1; r++){
						output->set( q*nMaxPxPerQuad + s*nPxPer2x1 + c*nRowsPer2x1 + r,
							copy->get( superrow+r, supercol+c ) );
					}
				}
			}
		}
		return 0;
	}


	//------------------------------------------------------------- createRawImageCSPAD
	// function to create 'raw' CSPAD images from plain 1D data, where
	// dim1 : 388 : rows of a 2x1
	// dim2 : 185 : columns of a 2x1
	// dim3 :   8 : 2x1 sections in a quadrant (align as super-columns)
	// dim4 :   4 : quadrants (align as super-rows)
	//
	//    +--+--+--+--+--+--+--+--+
	// q0 |  |  |  |  |  |  |  |  |
	//    |  |  |  |  |  |  |  |  |
	//    +--+--+--+--+--+--+--+--+
	// q1 |  |  |  |  |  |  |  |  |
	//    |  |  |  |  |  |  |  |  |
	//    +--+--+--+--+--+--+--+--+
	// q2 |  |  |  |  |  |  |  |  |
	//    |  |  |  |  |  |  |  |  |
	//    +--+--+--+--+--+--+--+--+
	// q3 |  |  |  |  |  |  |  |  |
	//    |  |  |  |  |  |  |  |  |
	//    +--+--+--+--+--+--+--+--+
	//     s0 s1 s2 s3 s4 s5 s6 s7
	//-------------------------------------------------------------
	int createRawImageCSPAD( const arraydata<double> *input, array2D<double> *&output ){
		if (!input){
			cerr << "Error in createRawImageCSPAD. No input. Aborting..." << endl;
			return 1;
		}

		delete output;
		output = new array2D<double>( nMaxQuads*nRowsPer2x1, nMax2x1sPerQuad*nColsPer2x1 );
		
		//sort into ordered 2D data
		for (int q = 0; q < nMaxQuads; q++){
			int superrow = q*nRowsPer2x1;
			for (int s = 0; s < nMax2x1sPerQuad; s++){
				int supercol = s*nColsPer2x1;
				for (int c = 0; c < nColsPer2x1; c++){
					for (int r = 0; r < nRowsPer2x1; r++){
						output->set( superrow+r, supercol+c, 
							input->get_atIndex( q*nMaxPxPerQuad + s*nPxPer2x1 + c*nRowsPer2x1 + r ) );
					}
				}
			}
		}
		
		//transpose to be conform with cheetah's convention
		output->transpose();
		return 0;
	}



	//------------------------------------------------------------- createAssembledImageCSPAD
	// expects 'pixX/Y' arrays to contain pixel count coordinates for each value in 'input'
	//-------------------------------------------------------------
	int createAssembledImageCSPAD( const arraydata<double> *input, const array1D<double> *pixX, const array1D<double> *pixY, array2D<double> *&output ){
		if (!input){
			cerr << "Error in createAssembledImageCSPAD. No input. Aborting..." << endl;
			return 1;
		}
		
		const double xmax = pixX->calcMax();
		const double xmin = pixX->calcMin();
		const double ymax = pixY->calcMax();
		const double ymin = pixY->calcMin();
		
		// calculate range of values in pixel arrays 
		// --> this will be the number of pixels in assembled image (add a safety margin)
		const int NX_CSPAD = (int) ceil(xmax-xmin) + 1;
		const int NY_CSPAD = (int) ceil(ymax-ymin) + 1;
		
		if (input->size() != pixX->size() || input->size() != pixY->size() ){
			cerr << "Error in array2D<double>::createAssembledImageCSPAD! Array sizes don't match. Aborting!" << endl;
			cerr << "size() of data:" << input->size() << ", pixX:" << pixX->size() << ", pixY:" << pixY->size() << endl;
			return 2;
		}else{
			//cout << "Assembling CSPAD image. Output (" << NX_CSPAD << ", " << NY_CSPAD << ")" << endl;
		}
		
		delete output;
		output = new array2D<double>( NY_CSPAD, NX_CSPAD );

		//sort into ordered 2D data
		for (int q = 0; q < nMaxQuads; q++){
			for (int s = 0; s < nMax2x1sPerQuad; s++){
				for (int c = 0; c < nColsPer2x1; c++){
					for (int r = 0; r < nRowsPer2x1; r++){
						int index = q*nMaxPxPerQuad + s*nPxPer2x1 + c*nRowsPer2x1 + r;

						//shift output arrays, if necessary, so that they start at zero
						output->set( (int)(pixY->get(index)-ymin), (int)(pixX->get(index)-xmin), input->get_atIndex(index) );
					}
				}
			}
		}
		
		return 0;
	}

}//namespace ns_cspad_util

