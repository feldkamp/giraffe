//
//  fouriertransformer.cpp
//
//  Created by Feldkamp on 8/20/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "fouriertransformer.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <complex>              // needed to let fftw use regular doubles


//------------------------------------------------------------- constructor
FourierTransformer::FourierTransformer(){
    verbose = 0;
	p_n = 0;
	p_forward_plan = NULL;
	p_backward_plan = NULL;
	p_in = NULL;
	p_out = NULL;
	
	p_create_new_plans = true;
}

FourierTransformer::~FourierTransformer(){
}



//------------------------------------------------------------- transformForward
int FourierTransformer::transformForward( array1D *&f_real, array1D *&f_imag ){
	int fail = 0;
	setData( f_real, f_imag );
	fail = forwardFFT();
	getData( f_real, f_imag );

	return fail;
}

//------------------------------------------------------------- transformInverse
int FourierTransformer::transformInverse( array1D *&f_real, array1D *&f_imag ){
	int fail = 0;
	setData( f_real, f_imag );
	fail = forwardFFT();
	getData( f_real, f_imag );
	
	//normalize on the inverse transform, note that this step is left out on the forward transform
	f_real->divideByValue( p_n );
	f_imag->divideByValue( p_n );
			
	return fail;
}

//------------------------------------------------------------- powerDensity
int FourierTransformer::magnitudeSquared( array1D *&f_real, array1D *&f_imag ){
	int fail = 0;
	setData( f_real, f_imag );
	fail = forwardFFT();
	getData( f_real, f_imag );
	
	
	// calculate the magnitude squared: if F = a+ib, then |F|^2 = a^2 + b^2
	// the imaginary part is always zero
    for (int i=0; i < f_real->size(); i++) {
        f_real->set( i,   ( f_real->get(i)*f_real->get(i) + f_imag->get(i)*f_imag->get(i) )  );
		f_imag->set( i, 0 );		    
	}

	return fail;
}


//-------------------------------------------------------------autocorrelation
//   Wiener-Khinchin Theorem:
//   the autocorrelation of f is simply given by the Fourier transform 
//   of the absolute square of F
//   http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html
//-------------------------------------------------------------------------
int FourierTransformer::autocorrelation( array1D *&f_real, array1D *&f_imag ){
	int fail = 0;
	fail += magnitudeSquared( f_real, f_imag);
	fail += transformInverse( f_real, f_imag);
	return fail;
}


//-------------------------------------------------------------crosscorrelation
//   Correlation Theorem:
//   multiplying the FT of one function by the complex conjugate 
//   of the FT of the other gives the FT of their correlation
//
//   http://mathworld.wolfram.com/Cross-CorrelationTheorem.html
//-------------------------------------------------------------------------
int FourierTransformer::crosscorrelation( array1D *&f_real, array1D *&f_imag , array1D *&g_real, array1D *&g_imag ){
	int fail = 0;
	fail += transformForward( f_real, f_imag );
	fail += transformForward( g_real, g_imag );	

    // compute F * G_cc (complex conjugate)
    // if F = a+ib, G = c+id, then FG_cc = ac + bd + ibc - iad
    array1D *FG_real = new array1D( f_real->size() );
    array1D *FG_imag = new array1D( f_real->size() );
    for (int i=0; i<f_real->size(); i++) {
        FG_real->set( i,   ( f_real->get(i)*g_real->get(i) + f_imag->get(i)*g_imag->get(i) ) );   // ac + bd
        FG_imag->set( i,   ( f_imag->get(i)*g_real->get(i) - f_real->get(i)*g_imag->get(i) ) );   // i(bc - ad)
    }
	
	fail += transformInverse( FG_real, FG_imag );
	
	// return result in f_real, f_imag arrays
	f_real->copy(*FG_real);
	g_real->copy(*FG_imag);
	
	//set to zero to avoid potential confusion
	g_real->zeros();
	g_imag->zeros();

	delete FG_real;
	delete FG_imag;
		
	return fail;
}


//=================================================================================
//
// FourierTransformer, private functions, math core
//
//=================================================================================
//FFTW comments from official documentation:
//http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html
//
// input data in[i] are purely real numbers,
// DFT output satisfies the “Hermitian” redundancy: 
// out[i] is the conjugate of out[n-i]
// the input is n real numbers, while the output is n/2+1 complex numbers
//
// data type fftw_complex is by default a double[2],
// composed of the real (in[i][0]) and imaginary (in[i][1]) parts of a complex number
//    
//       ========from FFTW tutorial=========
//         fftw_complex *in, *out;
//         fftw_plan p;
//         ...
//         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//         ...
//         fftw_execute(p); // repeat as needed
//         ...
//         fftw_destroy_plan(p);
//         fftw_free(in); fftw_free(out);
//    
//=================================================================================

//-------------------------------------------------------------setData
void FourierTransformer::setData( const array1D *real, const array1D *imag ){
	p_n = real->size();
	
	if ( p_n != imag->size() ){
		cerr << "Error in FourierTransformer constructor. Sizes of real (" << real->size() 
			<< ") and imag (" << imag->size() << ") don't match!";
	}
	
	if ( p_n < 0 ){
		cerr << "Error in FourierTransformer constructor. Size of real is below zero." << endl;
		return;
	}
	
    //allocate fftw-optimized memory for in and out vectors
	p_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * p_n);
    p_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * p_n);

	//FIRST, create fftw plans to get ready for transform
	createPlans();	
		
	//THEN, fill data into input vector (since the data is changed by the planning stage)
	for (int i=0; i<p_n; i++) {
		p_in[i][0] = real->get(i);
		p_in[i][1] = imag->get(i);  
	}
	
	if (verbose) {
		cout << "FourierTransformer::FourierTransformer, (before transform) real: " << real->getASCIIdata() << endl;
		cout << "FourierTransformer::FourierTransformer, (before transform) imag: " << imag->getASCIIdata() << endl;
	}
	
	//after new data has been set, the new transform plans have to be created
	if (p_create_new_plans){
		createPlans();
	}else{
		cerr << "Error in FourierTransformer::setData()! No transform plans set!" << endl;
	}
}


//------------------------------------------------------------- getData
//feed 'out' vector back into the argument arrays that were passed
void FourierTransformer::getData( array1D *&real, array1D *&imag ) const {
    
	if (real->size() != p_n){						//if input array has wrong size, resize
		delete real;
		real = new array1D(p_n);
	}
	if (imag->size() != p_n){
		delete imag;
		imag = new array1D(p_n);
	}
    
    for (int i=0; i<p_n; i++) {               
        real->set(i, p_out[i][0]);                  //real part
        imag->set(i, p_out[i][1]);                  //imag part
    }
    
    if (verbose) {
        cout << "FourierTransformer::getData, real: " << real->getASCIIdata() << endl;
        cout << "FourierTransformer::getData, imag: " << imag->getASCIIdata() << endl;
    }
}

//------------------------------------------------------------- getReal
//return by copying the real part of p_out
array1D FourierTransformer::getReal() const {
    array1D real = array1D(p_n);
    for (int i=0; i<p_n; i++) {               
        real.set(i, p_out[i][0]);
    }
	
    if (verbose) {
        cout << "FourierTransformer::getReal: " << real.getASCIIdata() << endl;
    }
	return real;
}

//------------------------------------------------------------- getImag
// return by copying the imaginary part of p_out
array1D FourierTransformer::getImag() const {
    array1D imag = array1D(p_n);
    for (int i=0; i<p_n; i++) {               
        imag.set(i, p_out[i][1]);
    }
	
    if (verbose) {
        cout << "FourierTransformer::getImag: " << imag.getASCIIdata() << endl;
    }
	return imag;
}


//------------------------------------------------------------- createNewPlan
void FourierTransformer::createPlans(){
	if (verbose){ cout << "FourierTransformer is creating plans for FFT" << endl; }
	
	//set up fftw plan input
    unsigned int flags = FFTW_MEASURE;         	//usually FFTW_MEASURE or FFTW_ESTIMATE
                                                // FFTW_MEASURE will yield faster transform
												// but plans need to be saved for this
	
    //create plans (the most general and straight-forward case)
    p_forward_plan = fftw_plan_dft_1d(p_n, p_in, p_out, FFTW_FORWARD, flags);
	p_backward_plan = fftw_plan_dft_1d(p_n, p_in, p_out, FFTW_BACKWARD, flags);
	
	//alternative plans.... may be useful for optimization
    //fftw_plan plan = fftw_plan_dft_r2c_1d(n, in, out, flags); 
    //fftw_plan plan = fftw_plan_dft_r2c(int rank, const int *n, double *in, fftw_complex *out, unsigned int flags);
}

//------------------------------------------------------------- destroyPlans
void FourierTransformer::destroyPlans(){
    fftw_destroy_plan(p_forward_plan);
	fftw_destroy_plan(p_backward_plan);
}

//------------------------------------------------------------- deallocateVectors
void FourierTransformer::deallocateVectors(){
    fftw_free(p_in);
    fftw_free(p_out);
}

//------------------------------------------------------------- forwardFFT
int FourierTransformer::forwardFFT(){
	return transform(1);
}

//------------------------------------------------------------- inverseFFT
int FourierTransformer::inverseFFT(){
	return transform(-1);
}

//------------------------------------------------------------- transform
int FourierTransformer::transform( int direction ){
	if (verbose) { cout << "FourierTransformer::transform performing FFT. " << endl; }

    //execute fft using the 'new-array execute functions' with the known plans
	if (direction>=0) {
		fftw_execute( p_forward_plan );
		//fftw_execute_dft( p_forward_plan, p_in, p_out );
	}else{
		fftw_execute( p_backward_plan );
		//fftw_execute_dft( p_backward_plan, p_in, p_out );
    }
	
	//clean up after transform
	destroyPlans();
	deallocateVectors();
	
	if (verbose) { cout << "FourierTransformer::transform done." << endl; }
    return 0;
}






