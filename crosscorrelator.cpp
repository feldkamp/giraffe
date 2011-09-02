/*
 *  crosscorrelator.cpp
 *
 *  XCCA analysis code 
 *  
 *  calculates the angular average of the intensity
 *  calculates the angular auto-correlation or cross-correlation
 *
 *  algorithm 1 (direct calculation) written by Jonas Sellberg, 2011-08-20 
 *  algorithm 2 (Fourier method) written by Jan Feldkamp, 2011-08-22
 *
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 */

#include "crosscorrelator.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cmath>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <cassert>

#include "fouriertransformer.h"
#include "arraydataIO.h"				// can be used, isn't be mandatory here


//=================================================================================
//
// constructors & destructor
//
//=================================================================================

//----------------------------------------------------------------------------constructor with C-style arrays
CrossCorrelator::CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, 
									int arraylength, int nphi, int nq1, 
									int nq2, int16_t *maskCArray){
    initPrivateVariables();
    
    //set basic properties for the size of the arrays
    setArraySize(arraylength);
	
	//check consistency of input arguments
	if (nphi > 0) {
		p_nPhi = nphi;
		p_nLag = (int) ceil(p_nPhi/2.0+1);		
	} else cerr << "ERROR in CrossCorrelator::constructor: nphi must be a positive integer." << endl;
	if (nq1 < 2) cerr << "ERROR in CrossCorrelator::constructor: nq1 must be a positive integer larger than 1." << endl;
	
	//check, user wants to run in 2D xaca or full 3D xcca mode
	if (nq2 == 0){
		setXccaEnable(false);
		p_nQ = nq1;
		p_autoCorrelation = new array2D( nQ(), nLag() );
	} else {
		setXccaEnable(true);
		p_nQ = (nq1 > nq2) ? nq1 : nq2;
		p_crossCorrelation = new array3D( nQ(), nQ(), nLag() );
	}
	
	//check, if a mask was given
	if (maskCArray == NULL){
		setMaskEnable( false );
	}else{
		setMaskEnable( true );
		setMask(maskCArray, arraylength);	
	}
	
    //copy data from array over to internal data structure
	setData(dataCArray, arraylength);

	//set q-calibration arrays (or generate a default)
	if ( qxCArray && qyCArray ){
		setQx(qxCArray, arraylength);
		setQy(qyCArray, arraylength);	
	}else{
		initDefaultQ();
	}

	initInternalArrays();
}


//----------------------------------------------------------------------------constructor with arraydata objects
CrossCorrelator::CrossCorrelator( array2D *dataArray, array2D *qxArray, array2D *qyArray, 
									int nphi, int nq1,
									int nq2, array2D *maskArray ){
	initPrivateVariables();
	
	if ( (dataArray->size() != qyArray->size()) || (dataArray->size() != qyArray->size()) ){
		cerr << "Warning in CrossCorrelator constructor! Array sizes don't match" << endl;
		cerr << "qxArray size = " << qxArray->size() << endl;
		cerr << "qyArray size = " << qyArray->size() << endl;
		cerr << "dataArray size = " << dataArray->size() << endl;
	}
	
    //set basic properties for the size of the arrays
    setArraySize(dataArray->size());
	p_nPhi = nphi;
	p_nLag = (int) ceil(p_nPhi/2.0+1);
	
	//check, user wants to run in 2D xaca or full 3D xcca mode
	if (nq2 == 0){
		setXccaEnable(false);
		p_nQ = nq1;
		p_autoCorrelation = new array2D( nQ(), nLag() );
	} else {
		setXccaEnable(true);
		p_nQ = (nq1 > nq2) ? nq1 : nq2;
		p_crossCorrelation = new array3D( nQ(), nQ(), nLag() );
	}

	//check, if a mask was given
	if (maskArray == NULL){
		setMaskEnable( false );
	}else{
		setMaskEnable( true );
		setMask(maskArray);	
	}
	
    //copy data from array over to internal data structure
    setData( dataArray );

	//set q-calibration arrays (or generate a default)
	if ( qxArray && qyArray ){
		setQx( qxArray );
		setQy( qyArray );	  
	}else{
		initDefaultQ();
	}
	
	initInternalArrays();
}


//----------------------------------------------------------------------------destructor
CrossCorrelator::~CrossCorrelator(){
	destroyInternalArrays();
	if (debug()){ 
		cout << "CrossCorrelator done after " << (clock()-p_creation_time)/(double)CLOCKS_PER_SEC << " seconds." << endl; 
	}
}


//=================================================================================
//
// initialize the class internals and set some defaults
//
//=================================================================================

//----------------------------------------------------------------------------initPrivateVariables
//make sure all private variables are initialized
//so that they don't contain or point to random memory
//----------------------------------------------------------------------------
void CrossCorrelator::initPrivateVariables(){
	p_arraySize = 1;
	p_qmin = 0;
	p_qmax = 0;
	p_deltaq = 0;
	p_phimin = 0;
	p_phimax = 0;
	p_deltaphi = 0;
	p_nQ = 0;
	p_nPhi = 0;
	p_nLag = 0;
	p_mask_enable = false;
	p_xcca_enable = true;
    p_outputdir = "";

	p_data = NULL;
	p_qx = NULL;
	p_qy = NULL;
	p_mask = NULL;
	p_polar = NULL;

	p_crossCorrelation = NULL;
	p_autoCorrelation = NULL;
	
	p_q = NULL;
	p_phi = NULL;	
	p_qAvg = NULL;
	p_iAvg = NULL;
	p_phiAvg = NULL;
	p_fluctuations = NULL;
	
	p_calculatePolarCoordinates = 0;
	p_calculateSAXS = 0;
	p_calculateXCCA = 0;
	p_updateDependentVariables = 0;
	
	p_mask_polar = NULL;
    p_table = NULL;
	
	p_qxmin = 0;
	p_qymin = 0;
	p_qxdelta = 0;
	p_qydelta = 0;
	
    p_debug = 0;
	p_creation_time = clock();
}

void CrossCorrelator::initInternalArrays(){
	//allocate all other internal objects
	p_table = new array2D(50, 50);
	
	p_q = new array1D(arraySize());
	p_phi = new array1D(arraySize());
	
	p_qAvg = new array1D(nQ());
	p_iAvg = new array1D(nQ());
	p_phiAvg = new array1D(nPhi());
}

void CrossCorrelator::destroyInternalArrays(){
	//free memory for objects
	delete p_data;
	delete p_qx;
	delete p_qy;
	delete p_mask;
    delete p_polar;
	
	delete p_mask_polar;
	delete p_table;
	
	delete p_q;
	delete p_phi;	
	delete p_qAvg;
	delete p_iAvg;
	delete p_phiAvg;
	delete p_fluctuations;
	
	delete p_crossCorrelation;
	delete p_autoCorrelation;
}


//----------------------------------------------------------------------------initDefaultQ
//if a q-calibration was not given to this class from the outside
//create a default one here
//----------------------------------------------------------------------------
void CrossCorrelator::initDefaultQ(){
    cout << "Initializing p_qx and p_qx vectors with default values." << endl;

    array1D *default_qx = new array1D(arraySize());
    array1D *default_qy = new array1D(arraySize());
    
    //set new values for deltaq and qmax
    p_deltaq = 1;                            
    p_qmax = arraySize()/2.*deltaq();

    for (int i=0; i<arraySize(); i++){
        default_qx->set(i, -qmax()+deltaq()*i );
        default_qy->set(i, -qmax()+deltaq()*i );
    }
    
    int onedim = (int) floor(sqrt(arraySize()));
    for (int i = 0; i < onedim; i++){
        for (int j = 0; j < onedim; j++){
            default_qx->set( i*onedim+j,    (i - (onedim-1)/2.)*deltaq() );
            default_qy->set( i*onedim+j,    (j - (onedim-1)/2.)*deltaq() );
        }
    }
	
	setQx(default_qx);
	setQy(default_qy);
		
	delete default_qx;
    delete default_qy;
}





//=================================================================================
//
// setters and getters for private variables
//
//=================================================================================

//----------------------------------------------------------------------------autoCorrelation
array2D *CrossCorrelator::autoCorr() const {
	return p_autoCorrelation;
}

//----------------------------------------------------------------------------crossCorrelation
array3D *CrossCorrelator::crossCorr() const {
	return p_crossCorrelation;
}


//----------------------------------------------------------------------------data
array1D *CrossCorrelator::data() const {
	return p_data;
}

void CrossCorrelator::setData( array1D *data ) {
	if (p_data) {
		delete p_data;
	}
	p_data = new array1D( *data );
	setArraySize( data->size() );
}

void CrossCorrelator::setData( array2D *data2D ) {
	if (p_data) {
		delete p_data;
	}
	p_data = new array1D( data2D );
	setArraySize( data2D->dim1() * data2D->dim2() );
}

void CrossCorrelator::setData( float *dataCArray, unsigned int size ){
	if (p_data) {
		delete p_data;
	}
	p_data = new array1D( dataCArray, size );
	setArraySize( size );
}

//----------------------------------------------------------------------------qx
array1D *CrossCorrelator::qx() const {
	return p_qx;
}

void CrossCorrelator::setQx( array1D *qx ){
	if (p_qx) {
		delete p_qx;
	}
	p_qx = new array1D( *qx );
}

void CrossCorrelator::setQx( array2D *qx ){
	if (p_qx) {
		delete p_qx;
	}
	p_qx = new array1D( qx );
}

void CrossCorrelator::setQx( float *qxArray, unsigned int size ){
	if (p_qx) {
		delete p_qx;
	}
	p_qx = new array1D( qxArray, size );
}


//----------------------------------------------------------------------------qy
array1D *CrossCorrelator::qy() const {
	return p_qy;
}

void CrossCorrelator::setQy( array1D *qy ) {
	if (p_qy) {
		delete p_qy;
	}
	p_qy = new array1D( *qy );
}

void CrossCorrelator::setQy( array2D *qy ) {
	if (p_qy) {
		delete p_qy;
	}
	p_qy = new array1D( qy );
}

void CrossCorrelator::setQy( float *qyArray, unsigned int size ){
	if (p_qy) {
		delete p_qy;
	}
	p_qy = new array1D( qyArray, size );
}

//----------------------------------------------------------------------------mask
array1D *CrossCorrelator::mask() const {
	return p_mask;
}

void CrossCorrelator::setMask( array1D *maskOneDArray ) {
	if (p_mask) {
		delete p_mask;
	}
	if (maskOneDArray) {
		p_mask = new array1D( *maskOneDArray );
		setMaskEnable( true );
		normalizeMask();
	}else{
		setMaskEnable( false );
	}
}

void CrossCorrelator::setMask( array2D *maskTwoDArray ) {
	if (p_mask) {
		delete p_mask;
	}
	if (maskTwoDArray) {
		p_mask = new array1D( maskTwoDArray );
		setMaskEnable( true );
		normalizeMask();
	}else{
		setMaskEnable( false );
	}
}

void CrossCorrelator::setMask(int16_t *maskCArray, unsigned int size){
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

void CrossCorrelator::normalizeMask(){
	//if the mask has a value other than zero (any value, pos. or neg.), set that to exactly 1
	for (int i = 0; i < mask()->size(); i++){
		double tolerance = 0.000001;	// float comparison, so allow for some error
		if ( mask()->get_atIndex(i) < tolerance && mask()->get_atIndex(i) > -tolerance ){
			mask()->set_atIndex(i, 0.);
		}else{
			mask()->set_atIndex(i, 1.);
		}
	}
}

void CrossCorrelator::setLookupTable( array2D *LUT ){
	calcLookupTableVariables( LUT->dim1(), LUT->dim2() );
	p_table->copy(*LUT);								//store a copy locally
}

void CrossCorrelator::setLookupTable( const int *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 ){
	calcLookupTableVariables( LUT_dim1, LUT_dim2 );	
	unsigned int tablelength = LUT_dim1*LUT_dim2;
	p_table->setDim1(LUT_dim1);						//set dimensions
	p_table->setDim2(LUT_dim2);
	p_table->arraydata::copy(cLUT, tablelength);		//store a copy in 'table' locally
}


array2D *CrossCorrelator::fluctuations() const {
	return p_fluctuations;
}


array2D *CrossCorrelator::polar() const {
	return p_polar;
}


array2D *CrossCorrelator::mask_polar() const {
	return p_mask_polar;
}

array2D *CrossCorrelator::lookupTable() const{
	return p_table;
}

//----------------------------------------------------------------------------arraySize
int CrossCorrelator::arraySize() const{
	return p_arraySize;
}

void CrossCorrelator::setArraySize( int arraySize_val ){
	if (arraySize_val > 0)
		p_arraySize = arraySize_val;
	else cerr << "ERROR in CrossCorrelator::setArraySize(): arraylength must be a positive integer." << endl;
}

//----------------------------------------------------------------------------matrixSize
int CrossCorrelator::matrixSize() const{
	return (int)floor(sqrt( (double)arraySize() ));            //assuming a square image
}

void CrossCorrelator::setMatrixSize( int matrixSize_val ){
	setArraySize( matrixSize_val*matrixSize_val );
}

//----------------------------------------------------------------------------qmin/qmax
void CrossCorrelator::setQminmax( double qmin_val, double qmax_val ){
	p_qmin = qmin_val;
	p_qmax = qmax_val;
	updateDependentVariables();
}

double CrossCorrelator::qmin() const{
	return p_qmin;
}

void CrossCorrelator::setQmin( double qmin_val ){
	p_qmin = qmin_val;
	updateDependentVariables();
}

double CrossCorrelator::qmax() const{
	return p_qmax;
}

void CrossCorrelator::setQmax( double qmax_val ){
	p_qmax = qmax_val;
	updateDependentVariables();
}

//----------------------------------------------------------------------------phimin/phimax
// CURRENTLY NOT USABLE WITH ALGORITHM 1, correlations of 360 degrees are always calculated
double CrossCorrelator::phimin() const{
	return p_phimin;
}

void CrossCorrelator::setPhimin( double phimin_val ){
	p_phimin = phimin_val;
}

double CrossCorrelator::phimax() const{
	return p_phimax;
}

void CrossCorrelator::setPhimax( double phimax_val ){
	p_phimax = phimax_val;
}

//------------------------------------------------------------------------------qmaxCArray
double CrossCorrelator::qmax2CArray( float *qxCArray, float *qyCArray, int arraylength ) {
	double qmax = 0;
	for (int i=0; i<arraylength; i++) {
		double qtemp = (double) sqrt(qxCArray[i]*qxCArray[i] + qyCArray[i]*qyCArray[i]);
		if (qtemp > qmax) qmax = qtemp;
	}
	return qmax;
}

double CrossCorrelator::qmax1CArray( float *qCArray, int arraylength ) {
	double qmax = 0;
	for (int i=0; i<arraylength; i++) {
		double qtemp; 
		if (qCArray[i] > 0) qtemp = (double) qCArray[i];
		else qtemp = (double) -qCArray[i];
		if (qtemp > qmax) qmax = qtemp;
	}
	return qmax;
}

//----------------------------------------------------------------------------maskEnable
bool CrossCorrelator::maskEnable() const {
    return p_mask_enable;
}

void CrossCorrelator::setMaskEnable( bool enable ){
	p_mask_enable = enable;
}

//----------------------------------------------------------------------------maskEnable
bool CrossCorrelator::xccaEnable() const {
    return p_xcca_enable;
}

void CrossCorrelator::setXccaEnable( bool enable ){
	p_xcca_enable = enable;
}

//----------------------------------------------------------------------------outputdir
string CrossCorrelator::outputdir(){
    return p_outputdir;
}

void CrossCorrelator::setOutputdir( std::string dir ){
	const char lastchar = dir.at( dir.size()-1 );
	if( lastchar != '/' ){		//if last character is not a slash, append one
		dir += '/';
	}
    p_outputdir = dir;
}

//----------------------------------------------------------------------------debug
int CrossCorrelator::debug() const {
    return p_debug;
}

void CrossCorrelator::setDebug( int debuglevel ){
    p_debug = debuglevel;
}

    
double CrossCorrelator::deltaq() const{						//getter only, dependent variable
	return p_deltaq;
}

double CrossCorrelator::deltaphi() const{						//getter only, dependent variable
	return p_deltaphi;
}

int CrossCorrelator::nQ() const{						//getter only, dependent variable
	return p_nQ;
}

int CrossCorrelator::nPhi() const{						//getter only, dependent variable
	return p_nPhi;
}

int CrossCorrelator::nLag() const{					//getter only, dependent variable
	return p_nLag;
}

void CrossCorrelator::updateDependentVariables(){		//update the values that depend on qmax and matrixSize
	// OLD BINNING
//	p_deltaq = 2*qmax()/(matrixSize()-1);
//	p_nQ = int(qmax()/p_deltaq+1+0.001);
//	p_deltaphi = 2*atan(1/(2*(p_nQ-1.0)));
//	p_nPhi = int(2*round(M_PI/p_deltaphi)); // make sure p_nPhi is even (exclude 2PI)
//	p_deltaphi = (double) 2.0*M_PI/(p_nPhi); // make sure deltaphi samples exactly an interval of 2PI
//	p_nLag = (int) round(p_nPhi/2.0+1);
	
	if (qmax() <= qmin()) {
		cout << "WARNING: Qmax <= Qmin, dependent variables not updated." << endl;
		return;
	}
	
	p_deltaq = (qmax()-qmin())/(nQ()-1);	// make sure deltaq samples start and stop
	p_deltaphi = (double) 2.0*M_PI/(p_nPhi);	// make sure deltaphi samples exactly an interval of 2PI
	
	if (debug() >= 1) {
		cout << "updated dependent variables:" << endl;
		cout << "qmin: " << qmin() << ", qmax: " << qmax() << ", p_deltaq: " << p_deltaq << ", p_nQ: " << p_nQ << ", p_deltaphi: " << p_deltaphi << ", p_nPhi: " << p_nPhi << ", p_nLag: " << p_nLag << endl;
	}
	
	p_updateDependentVariables++;
}

double CrossCorrelator::getQavg(unsigned index) const {
	return p_qAvg->get(index);
}

double CrossCorrelator::getPhiavg(unsigned index) const {
	return p_phiAvg->get(index);
}

double CrossCorrelator::getIavg(unsigned index) const {
	return p_iAvg->get(index);
}

double CrossCorrelator::getAutoCorrelation(unsigned index1, unsigned index2) const {
	return p_autoCorrelation->get(index1,index2);
}

double CrossCorrelator::getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const {
	return p_crossCorrelation->get(index1,index2,index3);
}




//=================================================================================
//
// CROSS-CORRELATION ALGORITHM 1
//
//=================================================================================

//--------------------------------------------------------calculatePolarCoordinates
// calculates polar coordinates for each pixel from cartesian coordinate system
void CrossCorrelator::calculatePolarCoordinates(double start_q, double stop_q) {
	
	// sanity limit check
	if (!p_qmax && !stop_q) {
		cerr << "ERROR in CrossCorrelator::calculatePolarCoordinates: Need to setQmax($VALUE) or specify Q-limits as arguments before running calculatePolarCoordinates()" << endl;
		return;
	} else if (start_q || stop_q) {
		p_qmin = start_q;
		p_qmax = stop_q;
		updateDependentVariables();
	}
	
	// function order check
	if (p_updateDependentVariables) {
		
		// calculate |q| for each pixel and bin lengths with correct resolution
		for (int i=0; i<arraySize(); i++) {
			p_q->set(i, round(sqrt( (qx()->get(i)*qx()->get(i))+(qy()->get(i)*qy()->get(i)) ) / deltaq()) * deltaq() );
		}		
		
		// calculate phi for each pixel and bin angles with correct resolution
		for (int i=0; i<arraySize(); i++) {
			double phii = 0;
			double qxi = qx()->get(i);
			double qyi = qy()->get(i);
			
			// setup UHP
			if (qxi == 0) { // make sure that the column with qx = 0 has angle 0 (q = 0 is assumed to have phi = 0)
				phii = 0;
			} else {
				phii = atan(qxi/qyi); // If qy = 0 and qx != 0, atan gives the correct result, but only for the UHP! Need to add PI for all LHP!
			}
			
			// correct LHP by adding PI
			if (qyi < 0) {
				phii += M_PI;
			}
			
			if (phii < -deltaphi()/2) { // make sure the binned angle is between 0 and 2PI-deltaphi()
				phii += 2*M_PI;
			}
			
			if (debug() >= 3) {
				if (phii < 0) cout << "phii: " << phii << ", nAngle(phii): " << round(phii/deltaphi()) << endl;
				else if (phii > (nPhi()-1)*deltaphi()) cout << "phii: " << phii << ", nAngle(phii): " << round(phii/deltaphi()) << endl;
			}
			
			p_phi->set( i, round(phii/deltaphi()) * deltaphi() );
		}
		
		// calculate vector of output |q|
		for (int i=0; i<nQ(); i++) {
			p_qAvg->set( i, qmin()+i*deltaq() );
		}

		// calculate vector of output angles
		for (int i=0; i<nPhi(); i++) {
			p_phiAvg->set( i, phimin()+i*deltaphi() );
		}
		
		p_calculatePolarCoordinates++;
	} else {
		cout << "WARNING: dependent variables must be updated before polar coordinates are calculated." << endl;
		updateDependentVariables();
		if (p_updateDependentVariables) calculatePolarCoordinates();
		else cerr << "ERROR in CrossCorrelator::calculatePolarCoordinates: dependent variables were not updated properly prior to use." << endl;
	}
}



//----------------------------------------------------------------------------calculateSAXS
// calculates the angular average of the intensity as a function of |q|
void CrossCorrelator::calculateSAXS(double start_q, double stop_q) {
	
	// sanity limit check
	if (!p_qmax && !stop_q) {
		cerr << "ERROR in CrossCorrelator::calculateSAXS: Need to setQmax($VALUE) or specify Q-limits as arguments before running calculateSAXS()" << endl;
		return;
	} else if (start_q || stop_q) {
		p_qmin = start_q;
		p_qmax = stop_q;
		updateDependentVariables();
	}
	
	// function order check
	if (p_updateDependentVariables) {
		
		// create counter array
		array1D *counter = new array1D( nQ() );
		
		// if calculateSAXS() has already been used, free and recreate p_qAvg, p_iAvg
		if (p_calculateSAXS) {
			delete p_qAvg;
			delete p_iAvg;
			p_qAvg = new array1D(nQ());
			p_iAvg = new array1D(nQ());
		}
		
		if (debug() >= 1) printf("calculating angular average of the intensity...\n");
		
		if (!p_calculatePolarCoordinates && !p_calculateSAXS) {
			// calculate |q| for each pixel and bin lengths with correct resolution
			for (int i=0; i<arraySize(); i++) {
				p_q->set(i, round(sqrt( (qx()->get(i)*qx()->get(i))+(qy()->get(i)*qy()->get(i)) ) / deltaq()) * deltaq() );
			}
		}
		
		if (!p_calculatePolarCoordinates || p_calculateSAXS) {
			// calculate vector of output |q|
			for (int i=0; i<nQ(); i++) {
				p_qAvg->set( i, qmin()+i*deltaq() );
			}			
		}
		
		// angular sum for each |q|
		for (int i=0; i<arraySize(); i++) {
			if (!maskEnable() || mask()->get(i)) {
				int qIndex = (int) round((p_q->get(i)-qmin())/deltaq()); // the index in qAvg[] that corresponds to q[i]
				if (p_q->get(i) <= qmax() && p_q->get(i) >= qmin()) {
					p_iAvg->set( qIndex, p_iAvg->get(qIndex)+data()->get(i) );
					counter->set( qIndex, counter->get(qIndex)+1 );
				} else if (debug() >= 3) printf("POINT EXCLUDED! q: %4.2f, qmin: %4.2f, qmax: %4.2f, nQ: %d, qIndex: %d\n", p_q->get(i), qmin(), qmax(), nQ(), qIndex );
			}
		}
		
		// normalize by number of pixels
		for (int i=0; i<nQ(); i++) {
			if (counter->get(i)) p_iAvg->set( i, p_iAvg->get(i)/counter->get(i) );
		}
		
		if (debug() >= 2) {
			printf("angular average of the intensity:\n");
			for (int i=0; i<nQ(); i++)
				cout << "Q: " << p_qAvg->get(i) << ",   \t# pixels: " << counter->get(i) << ",\tI: " << p_iAvg->get(i) << endl;
		}
		
		// free counter array
		delete counter;
		
		p_calculateSAXS++;
	} else {
		cout << "WARNING: dependent variables must be updated before SAXS intensity is calculated." << endl;
		updateDependentVariables();
		if (p_updateDependentVariables) calculatePolarCoordinates();
		else cerr << "ERROR in CrossCorrelator::calculateSAXS: dependent variables were not updated properly prior to use." << endl;
	}
}



//----------------------------------------------------------------------------calculateXCCA
void CrossCorrelator::calculateXCCA(double start_q, double stop_q) {
	
	// sanity limit check
	if (!p_qmax && !stop_q) {
		cerr << "ERROR in CrossCorrelator::calculateXCCA: Need to setQmax($VALUE) or specify Q-limits as arguments before running calculateXCCA()" << endl;
		return;
	} else if (start_q || stop_q) {
		calculatePolarCoordinates(start_q, stop_q);
	}
	
	// function order check
	if (p_calculatePolarCoordinates) {
		
		// create array of the intensities in polar coordinates with the correct binning
		delete p_polar;
		p_polar = new array2D( nQ(), nPhi() );
		
		// create array of the intensity fluctuations in polar coordinates with the correct binning
		delete p_fluctuations;
		p_fluctuations = new array2D( nQ(), nPhi() );
		
		// create local arrays over pixel counts in polar coordinates
		array2D *pixelCount = new array2D( nQ(), nPhi() );
		array2D *pixelBool = new array2D( nQ(), nPhi() );
		
		// if calculateXCCA() has already been used, free and recreate p_crossCorrelation, p_autoCorrelation
		if (p_calculateXCCA) {
			if (xccaEnable()) {
				delete p_crossCorrelation;
				p_crossCorrelation = new array3D( nQ(), nQ(), nLag() );
			} else {
				delete p_autoCorrelation;
				p_autoCorrelation = new array2D( nQ(), nLag() );
			}
		}
		
		if (debug() >= 1) printf("calculating polar arrays...\n");
		
		// calculate polar arrays
		for (int i=0; i<arraySize(); i++) {
			if (!maskEnable() || mask()->get(i)) {
				int qIndex = (int) round((p_q->get(i)-qmin())/deltaq()); // the index in qAvg[] that corresponds to q[i]
				int phiIndex = (int) round((p_phi->get(i)-phimin())/deltaphi()); // the index in phiAvg[] that corresponds to phi[i]
				if (qIndex >= 0 && qIndex < nQ() && phiIndex >= 0 && phiIndex < nPhi()) { // make sure qIndex and phiIndex are not out of array bounds
					polar()->set(qIndex, phiIndex, polar()->get(qIndex, phiIndex) + data()->get(i) );
					pixelCount->set(qIndex, phiIndex, pixelCount->get(qIndex,phiIndex)+1);
					if (pixelBool->get(qIndex, phiIndex) != 1) {
						pixelBool->set(qIndex, phiIndex, 1);
					}
				} else if (debug() >= 3) printf("POINT EXCLUDED! qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
			}
		}
		
		// calculate SAXS if not already calculated
		if (!p_calculateSAXS) calculateSAXS();
		
		// normalize polar array of the intensities and calculate fluctuations
		for (int i=0; i<nQ(); i++) {
			for (int j=0; j<nPhi(); j++) {
				if (pixelBool->get(i,j) != 0) {
					polar()->set(i, j, polar()->get(i,j)/pixelCount->get(i,j) );
					// subtract SAXS average for all shots or just subtract the SAXS from the specific shot?
					// second alternative is used here, since that always produces fluctuations with zero mean
					fluctuations()->set(i, j, polar()->get(i,j)-p_iAvg->get(i) );
				}
				if (debug() >= 2) printf("q: %f, phi: %f --> bool: %f, count: %f\n", p_qAvg->get(i), p_phiAvg->get(j), pixelBool->get(i, j), pixelCount->get(i, j));
			}
		}
		
		// calculate cross-correlation array and normalization array for cross-correlation	
		if (xccaEnable()) {
			
			if (debug() >= 1) printf("starting main loop to calculate cross-correlation...\n");
			
			for (int i=0; i<nQ(); i++) { // q1 index
				for (int j=0; j<nQ(); j++) { // q2 index
					for (int k=0; k<nLag(); k++) { // phi lag => phi2 index = (l+k)%nPhi
						double norm = 0;
						for (int l=0; l<nPhi(); l++) { // phi1 index
							int phi2_index = (l+k)%nPhi();
							crossCorr()->set(i,j,k, crossCorr()->get(i,j,k) + fluctuations()->get(i,l)*fluctuations()->get(j, phi2_index) );
							norm += pixelBool->get(i,l)*pixelBool->get(j, phi2_index);
						}
						if (norm) {
							crossCorr()->set(i, j, k, crossCorr()->get(i,j,k)/norm );
						} else { // fail code if no information exists about the specific element in the cross-correlation array
							crossCorr()->set(i, j, k, -2);
						}
					}
				} // for q2
			} // for q1
			
			// normalization loop for variances (since the diagonal elements w.r.t. Q of the cross-correlation array aren't calculated first)
			for (int i=0; i<nQ(); i++) { // q1 index
				for (int j=0; j<nQ(); j++) { // q2 index
					double variance1 = crossCorr()->get(i, i, 0);
					double variance2 = crossCorr()->get(j, j, 0);
					for (int k=0; k<nLag(); k++) { // phi lag => phi2 index = (l+k)%nPhi
//						if (p_iAvg->get(i) && p_iAvg->get(j)) {
//							// normalize the cross-correlation array with the SAXS intensity
//							crossCorr()->set(i, j, k, crossCorr()->get(i,j,k) / (p_iAvg->get(i)*p_iAvg->get(j)) );
//						} else { // fail code if average intensity is 0
//							crossCorr()->set(i, j, k, -1.5);
//						}
						if (variance1 && variance2) {
							// normalize by standard deviations (or the square root of the diagonal elements of the cross-correlation)
							crossCorr()->set(i, j, k, crossCorr()->get(i,j,k) / (sqrt(variance1)*sqrt(variance2)) );
						} else { // fail code if variances are 0
							crossCorr()->set(i, j, k, -1.5);
						}
					}
				} // for q2
			} // for q1
		} else {
			
			if (debug() >= 1) printf("starting main loop to calculate auto-correlation...\n");
			
			for (int i=0; i<nQ(); i++) { // q index
				double variance = 0;
				for (int k=0; k<nLag(); k++) { // phi lag => phi2 index = (l+k)%nPhi
					double norm = 0;
					for (int l=0; l<nPhi(); l++) { // phi1 index
						int phi2_index = (l+k)%nPhi();
						autoCorr()->set(i,k, autoCorr()->get(i,k) + fluctuations()->get(i,l)*fluctuations()->get(i, phi2_index) );
						norm += pixelBool->get(i,l)*pixelBool->get(i, phi2_index);
					}
					if (norm != 0) {
//						if (p_iAvg->get(i) != 0) {
//							// normalize the cross-correlation array with the SAXS intensity and the calculated normalization constant					
//							autoCorr()->set(i, k, autoCorr()->get(i,k) / (norm*p_iAvg->get(i)*p_iAvg->get(i)) );
//						} else { // fail code if average intensity is 0
//							autoCorr()->set(i, k, -1.5);						
//						}
						if (k == 0) {
							variance = autoCorr()->get(i, 0)/norm;
						}
						if (variance != 0) {
							// normalize by variance (or the zeroth element of the correlation)
							autoCorr()->set(i, k, autoCorr()->get(i,k) / (norm*variance) );
						} else { // fail code if variance is 0
							autoCorr()->set(i, k, -1.5);
						}
					} else { // fail code if no information exists about the specific element in the cross-correlation array
						autoCorr()->set(i, k, -2);
					}
				}		
			} // for q
		}
		
		delete pixelBool;
		delete pixelCount;
		
		if (debug() >= 1) printf("done calculating cross-correlation...\n");
		
		p_calculateXCCA++;
	} else {
		cout << "WARNING: polar coordinates must be calculated before XCCA is calculated." << endl;
		calculatePolarCoordinates();
		if (p_calculateXCCA) calculatePolarCoordinates();
		else cerr << "ERROR in CrossCorrelator::calculateXCCA: polar coordinates were not calculated properly prior to use." << endl;
	}	
}




//=================================================================================
//
// CROSS-CORRELATION ALGORITHM 2 (_FAST)
//
//=================================================================================


//----------------------------------------------------------------------------createLookupTable
// re-arrange data into a fast lookup table to get values fast using 
// val=lookup(x,y)
// dimensions of the argument 'table' determines the accuracy of the lookup
//----------------------------------------------------------------------------
int CrossCorrelator::createLookupTable( int lutNy, int lutNx ){
	if (debug()){ cout << "CrossCorrelator::createLookupTable() begin." << endl; }
    int retval = 0;
	
	calcLookupTableVariables( lutNy, lutNx );
	array2D *myTable = new array2D( lutNy, lutNx );
	
	//initialize table with -1, which can never be a real index
	const double initval = -1;
	myTable->addValue(initval);			
    
    if ( qx()->size()!=qy()->size() ) {
        cerr << "Error in createLookupTable! Array sizes don't match: " 
            << "qx=" << qx()->size() << ", qy=" << qy()->size() << endl;
        retval++;
    } else {
        double ix = 0;
        double iy = 0;
        for (int i = 0; i < qx()->size(); i++){           //go through all the data
            //get q-values from qx and qy arrays
            //and determine at what index (ix, iy) to put them in the lookup table
			double qxi = qx()->get(i);
			double qyi = qy()->get(i);
            ix = (qxi-p_qxmin) / p_qxdelta;
            iy = (qyi-p_qymin) / p_qydelta;
            
            //and fill table at the found coordinates with the data index
            //overwriting whatever value it had before
		    //(the add-one-half->floor trick is to achieve reasonable rounded integers)
            myTable->set( (int) floor(iy+0.5), (int) floor(ix+0.5), (double) i );
            
            /////////////////////////////////////////////////////////////////////////////////////////
            //ATTENTION: THIS METHOD WILL LEAD TO A LOSS OF DATA,
            //ESPECIALLY FOR SMALL TABLE SIZES,
            //BUT IT WILL BUY A LOT OF SPEED IN THE LOOKUP PROCESS
            //--> this should be improved to a more precise version, 
            //    maybe even one that allows the lookup(x,y) to interpolate
            //    for that to work, we need to find the four closest points in the data or so
            //    (for instance, instead of one index, the table could contain 
            //    a vector of all applicable indices)
            /////////////////////////////////////////////////////////////////////////////////////////
        }//for
		
		//after all is done, set the zero element to zero 
		//to make sure a failure of lookup() doesn't result in an acutal value
		myTable->set(0, 0, 0);		
    }//if

	if (debug()>1){ cout << "table = " << myTable->getASCIIdata() << endl; }
    if (debug()){ cout << "CrossCorrelator::createLookupTable() done." << endl; }
	
	//sanity check: how many, if any places are still at the default value?
	int initcount = 0;
	for (int i = 0; i < myTable->size(); i++){
		if ( myTable->get_atIndex(i) == initval ){
			initcount++;
		}
	}
	if (initcount > 0) {
		cout << "==============================================================================================" << endl;
		cout << "Warning in createLookupTable! There are still " << initcount << " unassigned values (" 
			<< ((double)initcount)/myTable->size()*100 << "%) in the lookup table." << endl;
		cout << "==============================================================================================" << endl;
	}
	
	//replace current table
	setLookupTable(myTable);
	
	delete myTable;	
    return retval;
}



//----------------------------------------------------------------------------calcLookupTableVariables
// calculate some variables needed to fill the lookup table
// and update the private variables
// which are important parameters for the lookup() function
//----------------------------------------------------------------------------
void CrossCorrelator::calcLookupTableVariables( int lutNy, int lutNx ){
	p_qxmin = qx()->calcMin();
    p_qxmax = qx()->calcMax();
    double qx_range = fabs(p_qxmax - p_qxmin);
    p_qxdelta = qx_range/(double)(lutNx-1);
    
    p_qymin = qy()->calcMin();
    p_qymax = qy()->calcMax();
    double qy_range = fabs(p_qymax - p_qymin);    
    p_qydelta = qy_range/(double)(lutNy-1);    
    
	if (debug()>1) {
		cout << "qx: min=" << p_qxmin << ", max=" << p_qxmax << ", range=" << qx_range << ", p_qxdelta=" << p_qxdelta << endl;  
		cout << "qy: min=" << p_qymin << ", max=" << p_qymax << ", range=" << qy_range << ", p_qxdelta=" << p_qydelta << endl;                    
	}
}



//---------------------------------------------------------------------------- subtractSAXSmean
int CrossCorrelator::subtractSAXSmean(){
	for(int q_ct=0; q_ct < polar()->dim1(); q_ct++){
		
		//get one row (fixed q) out of the polar coordinates
		array1D *row = new array1D( polar()->dim1() );
		polar()->getRow( q_ct, row );		
		
		//calculate average intensity in that row
		double avg = 0;
		if (!maskEnable()){			//no bad pixels --> just take the average of the row
			avg = row->calcAvg();
		}else{							//mask --> leave out bad pixels
			double sum = 0.;
			double valid = 0;
			for (int i = 0; i < row->size(); i++) {
				sum += row->get_atIndex(i);
				valid += mask_polar()->get(i, q_ct);		// if mask has a 1 here, the point is valid
			}
			avg = (valid > 0) ? (sum/((double)valid)) : 0;
		}
	
		row->subtractValue(avg);
		polar()->setRow(q_ct, row);
		delete row;
	}
	return 0;
}

//----------------------------------------------------------------------------transform to polar coordinates
int CrossCorrelator::calculatePolarCoordinates_FAST(){
	return calculatePolarCoordinates_FAST( 0, nQ() );
}

int CrossCorrelator::calculatePolarCoordinates_FAST( double start_q, double stop_q ) {
	
	//in principle, these could be arbitrary, but the FFT approach currently assumes wrap-around data
	const int start_phi = 0;
	const int stop_phi = 360;
	
	//find smallest common q in all directions, use that as upper limit, if necessary
	double abs_q_max = HUGE_VALF;
	if ( abs_q_max > fabs(p_qxmax)) { abs_q_max = fabs(p_qxmax); }
	if ( abs_q_max > fabs(p_qymax)) { abs_q_max = fabs(p_qymax); }
	if ( abs_q_max > fabs(p_qxmin)) { abs_q_max = fabs(p_qxmin); }
	if ( abs_q_max > fabs(p_qymin)) { abs_q_max = fabs(p_qymin); }
    stop_q = abs_q_max<stop_q ? abs_q_max : stop_q;  	//take the smallest value of the two

	if( debug()>1 ){ 
		cout << "CrossCorrelator::calculatePolarCoordinates_FAST" << endl; 
		cout << "varying scattering vector q from " << start_q << " to " <<  stop_q << " in " << nQ() << " steps, "
			<< "and angle phi full circle in " << nPhi() << " steps." << endl;
	}
    int retval = 0;

	//create new array2D to store polar coordinate representation
	delete p_polar;
    p_polar = new array2D( nQ(), nPhi() );
	if (!p_polar){
		cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST. polar couldn't be allocated." << endl;
		return 1;
	}
	
	if ( maskEnable() ) {
		//create new array2D to store polar coordinate representation
		delete p_mask_polar;
		p_mask_polar = new array2D( nQ(), nPhi() );
		if (!p_mask_polar){
			cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST. p_mask_polar couldn't be allocated." << endl;
			return 1;
		}
		setMaskEnable(false);	//disable mask feature for the purpose of treating the mask itself...
		calculatePolarCoordinates_FAST( mask(), mask_polar(), nPhi(), start_phi, stop_phi, nQ(), start_q, stop_q );
		setMaskEnable(true);	// ...turn back on
	}
	
	int novalue_count = calculatePolarCoordinates_FAST( data(), polar(), nPhi(), start_phi, stop_phi, nQ(), start_q, stop_q );
	if ( novalue_count > 0 ){
		cout << "Couldn't assign a true value in " << novalue_count << " cases (" 
			<< ((double)novalue_count)/nQ()/nPhi()*100 << "%) in the polar coordinate image of " 
			<< nPhi() << " x " << nQ() << " pixels." << endl;
    }

//	subtractSAXSmean();
	return retval;
}


//----------------------------------------------------------------------------calculatePolarCoordinates_FAST
// returns the number of times the lookup has failed
int CrossCorrelator::calculatePolarCoordinates_FAST( array1D* cartesian1D, array2D* polar2D, 
													int number_phi, double start_phi, double stop_phi, 
													int number_q, double start_q, double stop_q) const {
	double xcoord = 0.;
    double ycoord = 0.;
    double data_value = 0.;
    double q = 0.;
    double p = 0.;
    int q_ct = 0;
    int p_ct = 0;
	int novalue_ct = 0;
	
	//leave out the upper boundary, e.g. 200 <= q < 400, 0 <= phi < 360
	double step_q = (stop_q - start_q)/((double)number_q-1);
    double step_phi = (stop_phi - start_phi)/((double)number_phi-1);
	
	if (step_q < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST -- step_r value " 
            << step_q << " is smaller than zero." << endl;
    if (step_phi < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST -- step_phi value " 
            << step_phi << " is smaller than zero." << endl;
    
	for(q = start_q, q_ct=0; q_ct < number_q; q+=step_q, q_ct++){                        // q: for all rings/q-values
        if(debug()){cout << "#" << q_ct << ",q=" << q << "  " << std::flush;}
		for(p = start_phi, p_ct=0; p_ct < number_phi; p+=step_phi, p_ct++){				// phi: go through all angles

            //find lookup coordinates
			xcoord = q * cos(p*M_PI/180.);
			ycoord = q * sin(p*M_PI/180.);

            //lookup that value in original scattering data
			try {
            	data_value = lookup( xcoord, ycoord, cartesian1D );
            } catch (int e) {
				//cout << "Exception: " << e << endl;
				data_value = 0;									//setting the 0 here will introduce false correlations!!!
				novalue_ct++;
				if ( e == 0){
					cerr << "Exception " << e << ": table not allocated." << endl;
				}else{					
					//if (e == -1){ cerr << "Exception " << e << ": either interpolation in lookup failed or missing data." << endl; }
					//cases e = 1-4: asked for (x,y) values that are larger or smaller than the actual image
					//fail silently... data value remains zero. 
					
					if ( maskEnable() ){		//include this point in the mask to ignore it later in the analysis
						mask_polar()->set( q_ct, p_ct, 1);
					}
				}
			}//catch
			
			//assign the new value
			polar2D->set( q_ct, p_ct, data_value);
		}//for p
		
		
		//if a mask is used, attempt linear interpolation
		if (maskEnable()){
			array1D *datarow = new array1D( polar2D->dim2() );
			polar2D->getRow(q_ct, datarow);
			array1D *maskrow = new array1D( mask_polar()->dim2() );
			mask_polar()->getRow(q_ct, maskrow);
			//cout << "datarow: " << datarow->getASCIIdata();
			//cout << "maskrow: " << maskrow->getASCIIdata();
			
			double maskSum = maskrow->calcSum();
			if (maskSum >= maskrow->size()){					// nothing is masked (mask is 1 everywhere) --> go to next q
				continue;
			}
			if (maskSum < 0.5){									// everything is masked --> set data to zero, go to next q
				datarow->zeros();
				polar2D->setRow(q_ct, datarow);
				continue;
			}
			
			vector<int> start_m;	//mask starts
			vector<int> stop_m;		//mask stops
			int prev_maskval = 1;	//initial value 1
			int p_max = maskrow->size();
			
			for(int p_ct = 0; p_ct < p_max; p_ct++ ){
				int maskval = (int) maskrow->get(p_ct);
				if (maskval == 0){								//inside masked out region
					if ( p_ct == 0 ){							//this row begins with mask --> this is a start
						start_m.push_back(p_ct);
					}else if (p_ct == p_max-1){					//this row ends with mask --> this is a stop
						stop_m.push_back(p_ct);
					}else if (prev_maskval != 0){				//just entered this masked region --> just missed a start
						start_m.push_back(p_ct-1);
					}
				}
				if (maskval == 1){								//outside masked out region
					if (prev_maskval != 1){						//just exited this region --> this is a stop
						stop_m.push_back(p_ct);
					}
				}
				prev_maskval = maskval;
			}//for p_ct
			
			// there should be an equal number of starts and stops now, at least one of each
			if (start_m.size() != stop_m.size() || start_m.empty() || stop_m.empty()){
				cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST!";
				cerr << " start_m.size()=" << start_m.size() << " != stop_m.size()=" << stop_m.size() << endl;
				throw;
			}
			
			//calculate interpolation
			for (int k = 0; k < start_m.size(); k++){		//for each pair of starts and stops
				int thisstart = start_m.at(k);
				int thisstop = stop_m.at(k);
				double slope = (datarow->get(thisstop) - datarow->get(thisstart)) / (double)(thisstop-thisstart);
				double offset = datarow->get(thisstart);
				
				for (int i = thisstart; i < thisstop; i++){
					double fill_val = offset + (i-thisstart)*slope;
					datarow->set(i, fill_val);
				}
				
			}//for each start/stop
			
			polar2D->setRow(q_ct, datarow);
			//cout << "datarow after interpolation: " << datarow->getASCIIdata() << endl;
			delete datarow;
			delete maskrow;
		}//if maskEnable
	}//for q
	
    if (debug()) {
		cout << "polarCoordinates done. dimensions=(" << polar()->dim1() << " x " << polar()->dim2() << ")" << endl;
	    //write output of the intermediate files? (all zero by default, turn on for debugging or whatever)
		if (false){ cout << "data: " << cartesian1D->getASCIIdata() << endl; }
		if (false){ cout << "polar: " << polar2D->getASCIIdata() << endl; }
		//if (false){ io->writeToTiff( outputdir()+"polar.tif", polar ); }      //(arraydataIO needed for this)      
	}
	
	return novalue_ct;		
}




//----------------------------------------------------------------------------lookup
double CrossCorrelator::lookup( double xcoord, double ycoord, array1D *dataArray ) const {

    double value = 0.;    //return data value at the given coordinates
    int index = 0;      //to do that, the index in the data is determined first

    double xc = (xcoord-p_qxmin) / ((double)p_qxdelta);
    double yc = (ycoord-p_qymin) / ((double)p_qydelta);

    int ix = (int) floor( xc + 0.5 );		// round to nearest integer
    int iy = (int) floor( yc + 0.5 );
	
	bool enable_warnings = false;
    if ( !lookupTable() ){
        cerr << "Error in lookup! No lookup table was allocated." << endl;
		throw 0;
    } else if ( ix < 0 ){
        if (enable_warnings){
			cerr << "Error in lookup! xcoord=" << xcoord << " is too small.";
        	cerr << "   ix=" << ix << " < 0)" << endl;
		}
		throw 1;
    } else if ( ix >= lookupTable()->dim2() ){
        if (enable_warnings){
	        cerr << "Error in lookup! xcoord=" << xcoord << " is too large.";
	        cerr << "   ix=" << ix << " >= " << lookupTable()->dim2() << " (table dimx)" << endl;
		}
		throw 2;
    } else if ( iy < 0 ){
        if (enable_warnings){
	        cerr << "Error in lookup! ycoord=" << ycoord << " is too small.";
    	    cerr << "   iy=" << iy << " < 0)" << endl;
		}
		throw 3;
    } else if ( iy >= lookupTable()->dim1() ){
        if (enable_warnings){
			cerr << "Error in lookup! ycoord=" << ycoord << " is too large.";
			cerr << "   iy=" << iy << " >= " << lookupTable()->dim1() << " (table dimy)" << endl;
    	}
		throw 4;
	} else {
	    //create lookup index from original coordinates (assuming the data is properly centered)
        index = (int) lookupTable()->get(iy, ix);
		if ( index > 0 ){
			value = dataArray->get( index );
		} else {
			// if the value returned was negative, that means there is no lookup value assigned to this pair of (ix, iy)
			// this could be (1) because the table is too large, (2) because there are gaps in the data, e.g. CSPAD
			// solution for (1): create a binned value from the surrounding pixels
			int valid = 0;			//number of valid indices
			double sum = 0.;		//sum of retrieved values (at valid indices)
			vector<int> indices;
			if (ix+1 < lookupTable()->dim2()) 
				indices.push_back( (int) lookupTable()->get(iy, ix+1) );
			if (ix-1 >= 0)
				indices.push_back( (int) lookupTable()->get(iy, ix-1) );
			if (iy+1 < lookupTable()->dim1()) 
				indices.push_back( (int) lookupTable()->get(iy+1, ix) );
			if (iy-1 >= 0)
				indices.push_back( (int) lookupTable()->get(iy-1, ix) );
			
			for (int i = 0; i < indices.size(); i++){
				int tmpindx = indices.at(i);
				if ( tmpindx > 0 ){
					valid++;
					sum += dataArray->get( tmpindx );
				}
			}
			if ( valid > 0 ){
				value = sum / ((double)valid);
			} else {
				if (debug()){
					cout << "WARNING in lookup()! Couldn't find a value in LUT for coordinates (" << xcoord << ", " << ycoord << "). " << endl;
				}
				// search region would have to be expanded in a complete implementation....... 
				// the way it is should be enough for the moment, though, if the table isn't way too big......
				throw -1;
			}
		}
		
		//-----------------------!!!DEBUG!!!
		//make a hot pixel out of the one that was just looked up
		//careful! this actually changes the data in the input array!!
		//---> therefore, uncomment one of the following line only for debugging
		//p_dataCArray[index] = 100000.;
		//p_data->set( index, 1000);
		//-----------------------!!!DEBUG!!!
    }

    if( debug()>2 ){
        cout << "lookup (" << xcoord << ", " << ycoord 
            << ") --> LUT: (xc,yc)=(" << xc << ", " << yc 
            << ") ==> (" << ix << ", " << iy << ") "
            << "--> index=" << index << ", --> val=" << value << endl;
    }
    return value;
}






//----------------------------------------------------------------------------calculateXCCA_FAST
//
//normalize, according to eq. (3) in Altarelli, et al., PRB, 82, 104207 (2010)
// general, eq. (3):  C = ( <I(q1, phi) * I(q2, phi+deltaphi)> - <I(q1)>*<I(q2)> ) / ( <I(q1)>*<I(q2)> )
// autocorr, eq. (1): C =    <I(q, phi) * I(q, phi+deltaphi)>  /  ( <I(q)>^2 )  -  1
//
int CrossCorrelator::calculateXCCA_FAST(){
	if(debug()>1){ 
		cout << "CrossCorrelator::calculateXCCA_FAST ("<< polar()->dim1() << ", " << polar()->dim2() << ")" << endl; 
	}
	
	//do some sanity checks first
	if (!polar()) { 
		cerr << "No polar coordinate representation found. Call calculatePolarCoordinates_FAST() first." << endl; 
		throw 1; 
	}
	
	if (polar()->dim1()==0 || polar()->dim2()==0) {
		cerr << "Polar coordinate array has ill dimensions! x: " << polar()->dim1() << ", y: " << polar()->dim2() << endl;
		throw 2; 
	}
	
	//correlate data in polar, write result to corr
	try {
		if (!xccaEnable()){											
			// autocorrelation only
			autocorrelateFFT( polar(), autoCorr() );
		}else{														
			// full-blown cross-correlationp
			cerr << "CROSS-CORRELATE 3D USING FFT. PROCEED WITH CAUTION! Some features may not be implemented or fully tested." << endl;
			crosscorrelateFFT( polar(), crossCorr() );
		}//if
	} catch(int e) {
		cerr << "Exception " << e << " thrown: ";
		switch (e) {
			case 1:
				cerr << "Single row 'f' could not be allocated.";
				break;
			case 2:
				cerr << "Could not get single row.";
				break;
			case 3:
				cerr << "Could not calculate correlation.";
				break;
			default:
				cerr << "Unkown exception.";
		}//switch
		cerr << " Aborting. " << endl;
		throw e;	//re-throw
	}

	if(debug()){ cout << endl << "CrossCorrelator::calculateXCCA_FAST done." << endl; }			
    return 0;
}




//----------------------------------------------------------------------------autocorrelateFFT_byRow
int CrossCorrelator::autocorrelateFFT(array2D *polar2D, array2D *corr2D) const {
	int retval, fail = 0;
		
	//calculate the auto-correlation for all rings
	for(int q_ct=0; q_ct < polar2D->dim1(); q_ct++){
			
		//get one row out of the polar coordinates (at fixed q)
		array1D *f = new array1D( polar2D->dim2() );
		if (!f){ throw 1; } 
		polar2D->getRow( q_ct, f );
		if (f->size() == 0) { throw 2; }
				
		if (debug()>1){ cout << "   #" << q_ct << ", f before FFT: " << f->getASCIIdata() << endl; }
		if (debug()){ cout << q_ct << " " << std::flush; }
		
		//perform autocorrelation --> compute via FFT
		// initialize imaginary part to zero --> f is real
		array1D *f_imag = new array1D(f->size());						

		FourierTransformer *ft = new FourierTransformer();
		fail = ft->autocorrelation( f, f_imag );
		
		delete ft;
		delete f_imag;
		if (fail){ cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform (forward) failed." << endl; throw 3; }
		
		//normalize to zero-correlation
		f->divideByValue( f->get(0) );
		
		//feed result into corr
		//setRowPart will set only the first half of the data in f, the second half is redundant
		corr2D->setRow( q_ct, f);
		
		if (debug()>1){ cout << "   #" << q_ct << ", f after FFT: " << f->getASCIIdata() << endl; }
		delete f;
		f = NULL;
	}//for
	return retval;
}



//----------------------------------------------------------------------------autocorrelateFFT_byRow
int CrossCorrelator::crosscorrelateFFT(array2D *polar2D, array3D *corr3D) const {
	int retval, fail = 0;

	//calculate the cross-correlation for all combination of rings
	for(int fq_ct=0; fq_ct < polar2D->dim1(); fq_ct++){
		
		//get one row out of the polar coordinates
		array1D *f = new array1D( polar2D->dim2() );
		if (!f){ throw 1; } 
		polar2D->getRow( fq_ct, f );
		if (f->size() == 0) { throw 2; }
		
		//inner loop!
		for(int gq_ct=0; gq_ct < polar2D->dim1(); gq_ct++){

			//get second row out of the polar coordinates
			array1D *g = new array1D( polar2D->dim2() );
			if (!f){ throw 1; } 
			polar2D->getRow( gq_ct, g );		
			if (g->size() == 0) { throw 2; }
					
			if (debug()>1){ 
				cout << "   #" << fq_ct << ", f before FFT: " << f->getASCIIdata() << endl; 
			}
			if (debug()){ cout << fq_ct << " " << std::flush; }
			
			//perform CROSS-correlation --> compute via FFT
			if (!f || !g) {
				cerr << "CrossCorrelator::crosscorrelateFFT. Input not allocated." << endl;
				return 1;
			}
			if (f->size()==0 || g->size()==0) {
				cerr << "CrossCorrelator::crosscorrelateFFT. Input has size zero." << endl;
				return 2;
			}
			
			array1D *f_real = new array1D( *f );
			array1D *f_imag = new array1D(f->size());				// initialized to zeros --> f is real
			array1D *g_real = new array1D( *g );
			array1D *g_imag = new array1D(g->size());				// initialized to zeros --> g is real   
			
			// calculate the cross correlation, result is returned in f_real, f_imag
			FourierTransformer *ft = new FourierTransformer();
			fail = ft->crosscorrelation( f_real, f_imag, g_real, g_imag );
			delete ft;
			
			if (fail){ cerr << "Error in CrossCorrelator::crosscorrelateFFT. Transform (forward) failed." << endl; fail = 3; }
			
			// return result in original argument arrays (not really interested in imaginary part right now)
			f->copy( *f_real );
			//g->copy( *f_imag );
			
			delete f_real;
			delete f_imag;    
			delete g_real;
			delete g_imag;

			//normalize to zero-correlation
			f->divideByValue( f->get(0) );
			//g->divideByValue( f->get(0) );
	
			//feed result into corr
			//corr3D->setRow( fq_ct, gq_ct, f );		//(no such function exists, yet --> revisit later)
			for (int i = 0; i < f->size(); i++){
				corr3D->set( fq_ct, gq_ct, i, f->get(i) );
			}
			
			if (debug()>1){
				cout << "   #" << fq_ct << ", f before FFT: " << f->getASCIIdata() << endl; 
			}
			
			delete g;
			g = NULL;
		}//for: inner loop
		
		delete f;
		f = NULL;
	}//for: outer loop
	
	return retval;
}

