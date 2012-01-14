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
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, 
						int nphi, int nq1, int nq2 = 0);
	CrossCorrelator( arraydata<double> *dataArray, arraydata<double> *qxArray, arraydata<double> *qyArray, 
						int nphi, int nq1, int nq2 = 0);
	~CrossCorrelator();	

    void initPrivateVariables();
	void initInternalArrays();
	void destroyInternalArrays();
    void initDefaultQ();
	
	//---------------------------------------------main interface, calls different algorithms
	// master_algorithm: (1-4) call a combination of the partial algorithms below
	//		(1) calc coordinates with alg. 1, calc correlation with alg. 1
	//		(2) calc coordinates with alg. 2, calc correlation with alg. 2
	//		(3) calc coordinates with alg. 1, calc correlation with alg. 2
	//		(4) calc coordinates with alg. 2, calc correlation with alg. 1
	void run(int master_algorithm = 1);
	void run(double start_q, double stop_q, int master_algorithm = 1);
	
	// alg_coords: (1) full data with interpolation
	//             (2) only lookup certain values, FAST, but no interpolation
	// alg_corr:   (1) direct correlation calculation, works on data with gaps, slow
	//             (2) correlation via Fourier transform, fast, strong artifacts on data with gaps
	void run(double start_q, double stop_q, int alg_coords, int alg_corr);
	
	
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

	template <class T>
	void setLookupTable( T *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 ){
		calcLookupTableVariables( LUT_dim1, LUT_dim2 );	
		delete p_table;
		p_table = new array2D(cLUT, LUT_dim1, LUT_dim2);
		}
		
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
	
	int calculateSAXS_FAST();

	// calculate the auto- or cross-correlation function from p_polar and store it in p_corr
    int calculateXCCA_FAST();

	// "worker" functions for the one above:	
	// compute correlations in a 2D polar coordinate matrix
	// both functions do not change internal data in the class
	int autocorrelateFFT(array2D *polar2D, array2D *corr2D) const;
	int crosscorrelateFFT(array2D *polar2D, array3D *corr3D) const;
	
	
	//---------------------------------------------setters & getters

	//---------------------------------------------getters for variables
	int nQ() const;
	int nPhi() const;
	int nLag() const;
	double deltaq() const;
	double deltaphi() const;	
	
	//---------------------------------------------getters for calculated arrays
	array1D *qAvg() const;
	array1D *phiAvg() const;
	array2D *pixelCount() const;
	array1D *iAvg() const;
	array2D *fluctuations() const;	// intensity fluctuations in polar coordinates produced by calculatePolarCoordinates()
	array2D *polar() const;			// intensities in polar coordinates produced by calculatePolarCoordinates()/calculatePolarCoordinates_FAST()
	array2D *mask_polar() const;	// mask in polar coordinates produced by calculatePolarCoordinates_FAST()
	array2D *autoCorr() const;
	double getAutoCorrelation(unsigned index1, unsigned index2) const;
	array3D *crossCorr() const;	
	double getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const;
	
	//---------------------------------------------setters & getters for input data
	array1D *data() const;
	void setData( arraydata<double> *data );
	template <class T>
		void setData( T *dataCArray, unsigned int size ){
			if (p_data) {
				delete p_data;
			}
			p_data = new array1D( dataCArray, size );
			setArraySize( size );
		}
	
	array1D *qx() const;
	void setQx( arraydata<double> *qx );
	template <class T>
		void setQx( T *qxArray, unsigned int size ){
			if (p_qx) {
				delete p_qx;
			}
			p_qx = new array1D( qxArray, size );
		}

	array1D *qy() const;
	void setQy( arraydata<double> *qy );
	template <class T>
		void setQy( T *qyArray, unsigned int size ){
			if (p_qy) {
				delete p_qy;
			}
			p_qy = new array1D( qyArray, size );
		}

	array1D *mask() const;
	void setMask( arraydata<double> *mask );
	template <class T>
		void setMask(T *maskCArray, unsigned int size){
			if (p_mask) {
				delete p_mask;
			}
			if (maskCArray) {
				p_mask = new array1D( maskCArray, size );
				setMaskEnable( true );
				normalizeMask();
			}else{
				setMaskEnable( false );
			}
		}
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
	void setQminQmax( double qmin_val, double qmax_val );	
	
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
	array2D *p_pixelCount;		// pixel counts in polar coordinates
	
	array2D *p_fluctuations;	// intensity fluctuations in polar coordinates
    
	// function trackers
    int p_tracker_calculatePolarCoordinates;	// tracker for calculatePolarCoordinates()
	int p_tracker_calculateSAXS;				// tracker for calculateSAXS()
	int p_tracker_calculateXCCA;				// tracker for calculateXCCA()
	
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

