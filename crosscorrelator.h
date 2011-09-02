/*
 *  CrossCorrelator.h
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 04/27/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */
#ifndef _crosscorrelator_h
#define _crosscorrelator_h

#include "arrayclasses.h"


#include <string>

class CrossCorrelator {   

public:
	//---------------------------------------------constructors & destructor
	//arguments:
	//   data(C)Array: scattering data in a one-dimensional C-style array or a two-dimensional array2D object
	//   qx(C)Array: array of pixel x values of the same length
	//   qy(C)Array: array of pixel y values of the same length
	//   arraylength: length of all input arrays
	//   nphi: number of phi values, becomes length of the calculated correlation
	//   nq1: number of q values for which to calculate correlation
	//   nq2 (optional): causes computation of full 3D cross-correlation of nphi*nq1*nq2 values
	//   mask(C)Array: array that contains a mask of pixels to use or disregard
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, 
						int nphi, int nq1, int nq2 = 0, int16_t *maskCArray = NULL );
	CrossCorrelator( array2D *dataArray, array2D *qxArray, array2D *qyArray, 
						int nphi, int nq1, int nq2 = 0, array2D *maskArray = NULL  );
	~CrossCorrelator();
	

    void initPrivateVariables();
	void initInternalArrays();
	void destroyInternalArrays();
    void initDefaultQ();
	
	//---------------------------------------------calculations (Jonas's way)
	void calculatePolarCoordinates(double start_q = 0, double stop_q = 0);
	void calculateSAXS(double start_q = 0, double stop_q = 0);
	void calculateXCCA(double start_q = 0, double stop_q = 0);
	
    //---------------------------------------------alternative approach (Jan's way)
    // some of these functions have the byname _FAST to distinguish them from the ones above 
    // (for lack of a better name and in the hope that they may be fast. We'll see...)

    
	// a lookup table can be created within this object, if none was available externally
	// (for cheetah, performance is better if this is calculated once in advance 
	// and then handed to the CrossCorrelator for each shot)
	int createLookupTable( int Ny, int Nx );
	void calcLookupTableVariables( int lutNy, int lutNx );
	array2D *lookupTable() const;  
	void setLookupTable( array2D *LUT );
	void setLookupTable( const int *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 );
		
    // looks up the value closest to xcoord, ycoord in the data
    double lookup( double xcoord, double ycoord, array1D *dataArray ) const;
	
	//normalize the scattering data with the average for that specific q
	int subtractSAXSmean();

    // 'calculatePolarCoordinates' creates a 2D pattern in polar coordinates (r vs. phi)
    int calculatePolarCoordinates_FAST();
    int calculatePolarCoordinates_FAST(double start_q, double stop_q );
	
	// "worker" function for the one above
	int calculatePolarCoordinates_FAST( array1D* image, array2D* polar2D, 
										int number_phi,  double start_phi, double stop_phi, 
										int number_q, double start_q, double stop_q ) const;
	
	// calculate the auto- or cross-correlation function from p_polar and store it in p_corr
    int calculateXCCA_FAST();

	// "worker" functions for the one above:	
	// compute correlations in a 2D polar coordinate matrix
	// both functions do not change internal data in the class
	int autocorrelateFFT(array2D *polar2D, array2D *corr2D) const;
	int crosscorrelateFFT(array2D *polar2D, array3D *corr3D) const;
	
	
	//---------------------------------------------setters & getters

	//---------------------------------------------getters for dependent variables
	int nQ() const;
	int nPhi() const;
	int nLag() const;
	double deltaq() const;
	double deltaphi() const;	
	
	//---------------------------------------------getters for calculated arrays
	double getQavg(unsigned index) const;
	double getPhiavg(unsigned index) const;
	double getIavg(unsigned index) const;
	array2D *fluctuations() const;	// intensity fluctuations in polar coordinates produced by calculateXCCA()
	array2D *polar() const;			// intensities in polar coordinates produced by calculateXCCA()/calculatePolarCoordinates_FAST()
	array2D *mask_polar() const;	// mask in polar coordinates produced by calculatePolarCoordinates_FAST()
	array2D *autoCorr() const;
	double getAutoCorrelation(unsigned index1, unsigned index2) const;
	array3D *crossCorr() const;	
	double getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const;
	
	//---------------------------------------------setters & getters for input data
	array1D *data() const;
	void setData( array1D *data );
	void setData( array2D *data );
	void setData( float *dataCArray, unsigned int size );
	
	array1D *qx() const;
	void setQx( array1D *qx );
	void setQx( array2D *qx );
	void setQx( float *qxCArray, unsigned int size );

	array1D *qy() const;
	void setQy( array1D *qy );
	void setQy( array2D *qy );
	void setQy( float *qyCArray, unsigned int size );

	array1D *mask() const;
	void setMask( array1D *mask );
	void setMask( array2D *mask );	
	void setMask( int16_t *maskCArray, unsigned int size );
	void normalizeMask();

	
	int arraySize() const;                              // returns private variable p_arraySize
	void setArraySize( int arraySize_val );
	int matrixSize() const;                             // returns sqrt(p_arraySize)
	void setMatrixSize( int matrixSize_val );
                                                        // matrixSize is just a little helper for now
                                                        // to come up with a q-calibration
                                                        // need to change this soon...
	double qmax() const;
	void setQmax( double qmax_val );
	double qmin() const;
	void setQmin( double qmin_val );
	void setQminmax( double qmin_val, double qmax_val );	
	
	double phimin() const;
	void setPhimin( double phimin_val );
	double phimax() const;
	void setPhimax( double phimax_val );
	
	// jas: qmaxCArray() could easily be rewritten to use array1D qx, qy instead of CArrays if preferable
	double qmax2CArray( float *qxCArray, float *qyCArray, int arraylength ); // calculates qmax from 2 CArrays
    double qmax1CArray( float *qCArray, int arraylength ); // calculates qmax from 1 CArray
	
    bool maskEnable() const;	
	void setMaskEnable( bool enable );
	
    bool xccaEnable() const;
	void setXccaEnable( bool enable );
	
	// general output directory for the class
    std::string outputdir();
    void setOutputdir( std::string dir );

	// control the amount of (commandline) talking while running, default: 0
    int debug() const;
    void setDebug( int debuglevel );

	
	
private:
	//-------------------------------------------------required for both algorithms
	int p_arraySize;
	bool p_mask_enable;			// enables/disables masking of bad pixels
	bool p_xcca_enable;			// enables/disables cross-correlations (if disabled only autocorrelations are calculated)
    std::string p_outputdir;  	// the output directory if anything is dumped from withing this class (default is working dir)

	array1D *p_data;			// data storage
	array1D *p_qx;				// pixel x coordinate
	array1D *p_qy;				// pixel y coordinate
	array1D *p_mask;			// mask used to remove bad pixels
	array2D *p_polar;			// intensities in polar coordinates
	
	array2D *p_autoCorrelation;
	array3D *p_crossCorrelation;
	
	int p_debug;
	clock_t p_creation_time;
	
	//-------------------------------------------------required for alg1
	double p_qmin;				// start q for correlation calculations
	double p_qmax;				// stop q for correlation calculations
	double p_deltaq;			// step size (bin length) in q-direction
	double p_phimin;			// start phi for correlation calculations (currently disabled, 0 is always used)
	double p_phimax;			// stop phi for correlation calculations (currently disabled, 360 is always used)
	double p_deltaphi;			// step size (bin length) in phi direction
	int p_nQ;					// formerly samplingLength
	int p_nPhi;					// formerly samplingAngle
	int p_nLag;					// formerly samplingLag

	array1D *p_q;				// magnitude of q-vector (1st dimension in correlation) for each pixel
	array1D *p_phi;				// angle (2nd dimension in correlation) for each pixel
	
	array1D *p_qAvg;			// vector of output magnitudes of q-vector
	array1D *p_iAvg;			// vector of output average intensities for magnitudes of q-vector
	array1D *p_phiAvg;			// vector of output angles
	
	array2D *p_fluctuations;	// intensity fluctuations in polar coordinates
	
	void updateDependentVariables();
    
	// function trackers
    int p_calculatePolarCoordinates;	// tracker for calculatePolarCoordinates()
	int p_calculateSAXS;				// tracker for calculateSAXS()
	int p_calculateXCCA;				// tracker for calculateXCCA()
	int p_updateDependentVariables;		// tracker for updateDependentVariables()
	
	//-------------------------------------------------required for alg2
	array2D *p_mask_polar;
	array2D *p_table;					//lookup table

	double p_qxmax;
	double p_qymax;	
	double p_qxmin;
	double p_qymin;
	double p_qxdelta;
	double p_qydelta; 
};





#endif

