//
//  fouriertransformer.h
//	an easy-to-use wrapper for the Fourier transform calculated in the FFTW3 library
//
//  Created by Feldkamp on 8/20/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//


#include "arrayclasses.h"
#include <fftw3.h>


class FourierTransformer{
public:
    FourierTransformer();
    ~FourierTransformer();	
	
	//each of these functions is passed the real and imaginary part of a complex function f (f_real, f_imag)
	//result is returned in these same arguments
	int transformForward( array1D *&f_real, array1D *&f_imag );					
	int transformInverse( array1D *&f_real, array1D *&f_imag );					// includes 1/N normalization
	int magnitudeSquared( array1D *&f_real, array1D *&f_imag );
	int autocorrelation( array1D *&f_real, array1D *&f_imag );
	int crosscorrelation( array1D *&f_real, array1D *&f_imag, array1D *&g_real, array1D *&g_imag );
	
private:
    int verbose;
	fftw_complex *p_in;				//internal complex input array
	fftw_complex *p_out;			//internal complex output array
	fftw_plan p_forward_plan;
	fftw_plan p_backward_plan;
	int p_n;						//size of the input (and output) arrays
	
	bool p_create_new_plans;		//default: new plans are created every time
									//alternatively, the user could be enabled to set plans manually (not implemented currently)
	
	void setData( const array1D *real, const array1D *imag );
	void getData( array1D *&real, array1D *&imag ) const;		// return within passed arguments
	//after the transform, use these functions to ask for the transformed data
	array1D getReal() const;									// return by copy (may be slower)
	array1D getImag() const;									// return by copy (may be slower)
	
	void createPlans();
	void destroyPlans();
	void deallocateVectors();
	
    //wrappers for the FFTW discrete Fourier Transform
	int forwardFFT();					// == transform(1)
	int inverseFFT();					// == transform(-1)
	int transform( int direction=1 );
};
