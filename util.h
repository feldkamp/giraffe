#ifndef UTIL_H
#define UTIL_H
//
//  util.h
//  -- a collection of useful utility functions that are not specific to a certain module
//
//  Created by Feldkamp on 9/28/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "arrayclasses.h"

namespace ns_cspad_util{
	const int nRowsPerASIC = 194;
	const int nColsPerASIC = 185;
	const int nRowsPer2x1 = nRowsPerASIC * 2;								// 388;
	const int nColsPer2x1 = nColsPerASIC;									// 185;
	const int nMax2x1sPerQuad = 8;
	const int nMaxQuads = 4;
	const int nPxPer2x1 = nColsPer2x1 * nRowsPer2x1;						// 71780;
	const int nMaxPxPerQuad = nPxPer2x1 * nMax2x1sPerQuad;					// 574240;
	const int nMaxTotalPx = nMaxPxPerQuad * nMaxQuads; 						// 2296960;  //(2.3e6)
	
	//-------------functions for the CSPAD detector-------------
	int create1DFromRawImageCSPAD( const array2D *input, array1D *&output );
	int createRawImageCSPAD( const arraydata<double> *input, array2D *&output );
	int createAssembledImageCSPAD( const arraydata<double> *input, const array1D *pixX, const array1D *pixY, array2D *&output );
}	
	
namespace ns_histogram_util{	
	int getHistogramInBoundaries( const arraydata<double> *data, array1<int> *&hist, array1<double> *&bins, unsigned int nBins, double min, double max );
	int getHistogram( const arraydata<double> *data, array1<int> *&hist, array1<double> *&bins, unsigned int nBins = 20 );
	std::string getHistogramInBoundariesASCII( const arraydata<double> *data, unsigned int nBins, double min, double max );
	std::string getHistogramASCII( const arraydata<double> *data, unsigned int nBins = 20 );
}


#endif
