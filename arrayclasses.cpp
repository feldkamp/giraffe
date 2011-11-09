/*
 *  ArrayClasses.cpp
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 07/22/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */

#include "arrayclasses.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <iomanip>

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <fstream>
using std::ofstream;

#include <vector>
using std::vector;

#include <cmath>
#include <cassert>
#include <exception>




//=================================================================================
//
// CLASS IMPLEMENTATION OF arraydata
//
//=================================================================================
arraydata::arraydata( unsigned int size_val ){
    init();
	p_size = size_val;
	try{
		p_data = new double[size_val];
	}catch (std::exception& e){
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;
		cerr << "Standard exception: " << e.what() << endl;
	}
	zeros();					// set all elements to zero initially
}

arraydata::arraydata( const arraydata *array ){
	init();
	this->copy( array );
}

arraydata::arraydata( const array1D *array ){
	init();
	this->copy( array );
}

arraydata::arraydata( const array2D *array ){
	init();
	this->copy( array );
}

arraydata::arraydata( const array3D *array ){
	init();
	this->copy( array );
}

arraydata::arraydata( const array4D *array ){
	init();
	this->copy( array );
}

arraydata::arraydata( const vector<double> vec ){
	init();
	this->copy( vec );
}


//copy constructor
arraydata::arraydata( const arraydata &src ){                       
    init();
    copy( src );
}

//assignment operator
arraydata & arraydata::operator=(const arraydata & src){
    if ( this != &src ){
        this->destroy();
        init();
        copy( src );
    }
    return *this;
}

arraydata::~arraydata(){
	this->destroy();	
}

//--------------------------------------------------------------------helper functions
void arraydata::init(){
    p_size = 0;
    p_data = NULL;
}

void arraydata::destroy(){
    delete []p_data;
	p_data = NULL;
}

double *arraydata::data() const{
    return (double *)p_data;
}

//-----------------------------------------------------copy
void arraydata::copy( const arraydata *src ){										//src type: arraydata pointer
    this->copy( src->data(), src->size());
}

void arraydata::copy( const array1D *src ){											//src type: array1D pointer
	this->copy( src->data(), src->size() );
}

void arraydata::copy( const array2D *src ){											//src type: array2D pointer
	this->copy( src->data(), src->size() );
}

void arraydata::copy( const array3D *src ){											//src type: array3D pointer
	this->copy( src->data(), src->size() );
}

void arraydata::copy( const array4D *src ){											//src type: array4D pointer
	this->copy( src->data(), src->size() );
}
	
void arraydata::copy( const arraydata& src ){										//src type: arraydata object
    this->copy( src.data(), src.size());
}

void arraydata::copy( const std::vector<double> vec ){
	p_size = (unsigned int) vec.size();
	if (p_size > 0){
		this->destroy();
		p_data = new double[ p_size ];
		for (unsigned int i = 0; i < p_size; i++) {
			p_data[i] = (double)vec.at(i);
		}
	}else{
		p_data = NULL;
	}
}

//--------------------------------------------------------------------
//time-critical function, check with assert is disabled in release configuration
double arraydata::get_atIndex( unsigned int i) const{
	assert(p_data);
	assert(i < size());
	return p_data[i];
}


//--------------------------------------------------------------------
//time-critical function, check with assert is disabled in release configuration
void arraydata::set_atIndex( unsigned int i, double val){
	assert(i<size());
	p_data[i] = val;
}


//--------------------------------------------------------------------
unsigned int arraydata::size() const{
	return p_size;
}


//--------------------------------------------------------------------
void arraydata::zeros(){					//set all elements to zero
	for (int i = 0; i < size(); i++) {
		set_atIndex(i, 0);
	}
}

void arraydata::zeros( unsigned int start, unsigned int stop ){
	if ( start >= stop || stop > size() ){
		cerr << "Error in arraydata::zero("<< start << ", " << stop << "). Check boundaries." << endl;
		throw;
	}
	for (int i = start; i < stop; i++) {
		set_atIndex(i, 0);
	}
}

//--------------------------------------------------------------------
void arraydata::ones(){					//set all elements to one
	for (int i = 0; i < size(); i++) {
		set_atIndex(i, 1);
	}	
}

void arraydata::ones( unsigned int start, unsigned int stop ){
	if ( start >= stop || stop > size() ){
		cerr << "Error in arraydata::ones("<< start << ", " << stop << "). Check boundaries." << endl;
		throw;
	}
	for (int i = start; i < stop; i++) {
		set_atIndex(i, 1);
	}
}


void arraydata::range( double neg, double pos ){	//set elements to a range of values, given by the boundaries
	double delta = (pos-neg)/(size()-1);
	for (int i = 0; i < size(); i++) {
		set_atIndex(i, neg+delta*i);
	}	
}

//-------------------------------------------------------------------- calcMin
// ignores the position of the minimum
double arraydata::calcMin() const{
	int pos = 0;					
	return this->calcMin(pos);
}

//-------------------------------------------------------------------- calcMin
// returns the position of the minimum
double arraydata::calcMin(int &pos) const{
	double tempmin = +INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atIndex(i) < tempmin) {
			tempmin = get_atIndex(i);
			pos = i;
		}
	}
	return tempmin;
}

//-------------------------------------------------------------------- calcMax
// ignores the position of the minimum
double arraydata::calcMax() const{
	int pos = 0;					
	return this->calcMax(pos);
}

//-------------------------------------------------------------------- calcMax
// returns the position of the maximum
double arraydata::calcMax(int &pos) const{
	double tempmax = -INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atIndex(i) > tempmax) {
			tempmax = get_atIndex(i);
			pos = i;
		}
	}
	return tempmax;
}

//------------------------------------------------------------- calcSum
double arraydata::calcSum() const{
	double sum = 0.;
	for (int i = 0; i < size(); i++) {
		sum += get_atIndex(i);
	}	
	return sum;
}

//------------------------------------------------------------- calcAvg
double arraydata::calcAvg() const{
	double avg = this->calcSum() / ((double)size());
	return avg;
}

//------------------------------------------------------------- getASCIIdataAsRow
string arraydata::getASCIIdataAsRow() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "Array data has size zero." << endl;
	}else{
		for (int i = 0; i<size(); i++) {
			osst << get_atIndex(i) << " ";
		}
		osst << endl;
	}
    return osst.str();
}

//------------------------------------------------------------- getASCIIdataAsColumn
string arraydata::getASCIIdataAsColumn() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "Array data has size zero." << endl;
	}else{
		for (int i = 0; i<size(); i++) {
			osst << get_atIndex(i) << endl;
		}
	}
    return osst.str();
}


//--------------------------------------------------------------------array math
//multiply each element by a numerical factor
int arraydata::addValue( double val ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i) + val);
    }
    return 0;
}

//multiply each element by a numerical factor
int arraydata::subtractValue( double val ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i) - val);
    }
    return 0;
}

//multiply each element by a numerical value
int arraydata::multiplyByValue( double value ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i) * value);
    }
    return 0;
}

//divide each element by a numerical value (enforcing float division)
int arraydata::divideByValue( double value ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, ((double) this->get_atIndex(i)) / value);
    }
    return 0;
}


//add each element by an element from a second array
int arraydata::addArrayElementwise( const arraydata *secondArray ){
    if (!secondArray){
		cerr << "Error in arraydata::addArrayElementwise! Second array not allocated. ";
		return 2;
	}
	if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::addArrayElementwise! Array sizes don't match. ";
        cerr << "(this array size " << this->size() << " != second array size " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i)+secondArray->get_atIndex(i));
    }
    return 0;
}

//subtract each element by an element from a second array
int arraydata::subtractArrayElementwise( const arraydata *secondArray ){
    if (!secondArray){
		cerr << "Error in arraydata::subtractArrayElementwise! Second array not allocated. ";
		return 2;
	}
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::subtractArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i)-secondArray->get_atIndex(i));
    }
    return 0;
}

//multiply each element by an element from a second array
int arraydata::multiplyByArrayElementwise( const arraydata *secondArray ){
    if (!secondArray){
		cerr << "Error in arraydata::multiplyArrayElementwise! Second array not allocated. ";
		return 2;
	}
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::multiplyArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i) * secondArray->get_atIndex(i) );
    }
    return 0;
}

//divide each element by an element from a second array
int arraydata::divideByArrayElementwise( const arraydata *secondArray ){
    if (!secondArray){
		cerr << "Error in arraydata::divideArrayElementwise! Second array not allocated. ";
		return 2;
	}
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::divideArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atIndex(i, this->get_atIndex(i)/secondArray->get_atIndex(i));
    }
    return 0;
}


//------------------------------------------------------------- applyMask
int arraydata::applyMask( arraydata* mask, double checkval, double rejectval ){
	if (size() != mask->size()){
		cerr << "Error in arraydata::applyMask! Data size " << size() << " doesn't match mask size " << mask->size() << "." << endl;
		return 1;
	}
	for (int i = 0; i < size(); i++){
		if (mask->get_atIndex(i) < checkval){
			//replace by fill value
			set_atIndex(i, rejectval);
		}else{
			//do nothing, keep the original data
		}
	}	
}

//------------------------------------------------------------- histogram
int arraydata::getHistogramInBoundaries( array1D *&hist, array1D *&bins, unsigned int nBins, double min, double max ){
	double range = max-min;
	double binwidth = range/(double)(nBins);
	double bindelta = range/(double)(nBins-1);
	//cout << "min: " << min << ", max:" << max << ", range: " << range << ", binwidth: " << binwidth << ", nBins: " << nBins << endl;
	
	if (range == 0 && hist->size() > 1){
		cerr << "WARNING in arraydata::getHistogramInBoundaries. max == min == " << max << endl;
	}
	
	//make bins array
	delete bins;
	bins = new array1D(nBins);
	for (int b = 0; b < bins->size(); b++){
		bins->set( b, min + b*binwidth );
	}
	
	//make histogram
	delete hist;
	hist = new array1D(nBins);
	for (int i = 0; i < this->size(); i++){
		unsigned int binnum = (unsigned int) floor( (this->get_atIndex(i)-min) / (double)bindelta );
		if (binnum < hist->size()){
			hist->set( binnum, hist->get(binnum) + 1 );
		}
	}
	return 0;
}


std::string arraydata::getHistogramInBoundariesASCII( unsigned int nBins, double min, double max ){
	array1D *hist = new array1D();
	array1D *bins = new array1D();
	
	this->getHistogramInBoundaries( hist, bins, nBins, min, max );

	int maxpos = 0;
	double histmax = hist->calcMax( maxpos );
	int nMarkers = 20;
	int numberPerMarker = (int) floor(histmax/nMarkers);
	if (numberPerMarker < 1){
		numberPerMarker = 1;
	}
	
	double binsize = (max-min)/(double)(nBins+1);
	ostringstream osst;
	osst << std::setw(8) << "#" << " (" << std::setw(12) << "low" << ", " 
		<< std::setw(12) << "high" << "): " << std::setw(12) << "occurrence" << endl;
	osst << "---------------------------------------------------" << endl;
	for (int i = 0; i < bins->size(); i++){
		osst << std::setw(8) << i << " (" << std::setw(12) << bins->get(i) << ", " 
			<< std::setw(12) << bins->get(i)+binsize << "): " 
			<< std::setw(12) << hist->get(i) << "  |";
		for (int m = 0; m < (hist->get(i)/(double)numberPerMarker); m++){ osst << "*"; }
		osst << endl;
	}
	osst << "---------------------------------------------------" << endl;
	osst << "max:" << std::setw(4) << maxpos << " (" << std::setw(12) << bins->get(maxpos) << ", " 
		<< std::setw(12) << bins->get(maxpos)+binsize << "): " << std::setw(12) << hist->get(maxpos) << "  |" << endl;
			
	delete hist;
	delete bins;
	return osst.str();
}

int arraydata::getHistogram( array1D *&hist, array1D *&bins, unsigned int nBin ){
	double min = this->calcMin();
	double max = this->calcMax();
	return this->getHistogramInBoundaries( hist, bins, nBin, min, max );
}

std::string arraydata::getHistogramASCII( unsigned int nBin ){
	double min = this->calcMin();
	double max = this->calcMax();
	return this->getHistogramInBoundariesASCII( nBin, min, max );
}





//=================================================================================
//
// CLASS IMPLEMENTATION OF array1D
//
//=================================================================================
//-----------------------------------------------------constructors & destructors
array1D::array1D( unsigned int size_dim1 ) 
		: arraydata(size_dim1){
    setDim1( size_dim1 );
}

array1D::array1D( arraydata *array )
		: arraydata( array ) {
    setDim1( array->size() );
}

array1D::array1D( array1D *arrayOneD )
		: arraydata( arrayOneD ) {
    setDim1( arrayOneD->size() );
}

//constructor to generate a 1D array from a 2D array
array1D::array1D( array2D* arrayTwoD ) 
        : arraydata( arrayTwoD->size() ){
 	setDim1( arrayTwoD->size() );

    //copy contents of dataTwoD to this array1D object
    p_size = arrayTwoD->size();
    if (p_size > 0){
        for (int i = 0; i < p_size; i++) {
            p_data[i] = arrayTwoD->get_atIndex(i);
        }
    }else{
        p_data = NULL;
    }
}


//-----------------------------------------------------destructor
array1D::~array1D(){
}


//-----------------------------------------------------copy
void array1D::copy( const array1D& src ){
    setDim1( src.size() );
    this->arraydata::copy( src.data(), src.size() );
}



//-----------------------------------------------------get
//time-critical function, check with assert is disabled in release configuration
double array1D::get( unsigned int i ) const{
	assert(i < dim1());
	return arraydata::get_atIndex(i);		
}

//-----------------------------------------------------set
//time-critical function, check with assert is disabled in release configuration
void array1D::set( unsigned int i, double value ){
	assert(i < dim1());
	arraydata::set_atIndex(i, value);
}


//-----------------------------------------------------setters & getters
unsigned int array1D::dim1() const{
	return p_dim1;
}

void array1D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}


//------------------------------------------------------------- getASCIIdata
string array1D::getASCIIdata( bool annotate ) const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "1D data has size zero." << endl;
	}else{
		if (annotate){
			osst << "1D data, dim1=" << dim1() << ", size=" << size() << endl;
			osst << " [";
		}
		for (int i = 0; i<dim1(); i++) {
			osst << " " << get(i);
		}
		if (annotate){osst << "]";}
		osst << endl;
	}
    return osst.str();
}











//=================================================================================
//
// CLASS IMPLEMENTATION OF array2D
//
//=================================================================================
//-----------------------------------------------------constructors & destructors
array2D::array2D( unsigned int size_dim1, unsigned int size_dim2 )
		: arraydata(size_dim1*size_dim2){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
}

array2D::array2D( arraydata* data, unsigned int size_dim1, unsigned int size_dim2 )
        : arraydata( data ){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
    if (data->size() != size_dim1*size_dim2) {
        cerr << "WARNING in array2D::array2D. Inconsistent array size. ";
        cerr << "size1D=" << data->size() << ", size2D=" << size_dim1*size_dim2 
			<< "=" << size_dim1 << "*" << size_dim2 << "" << endl;
    }
}

//constructor to generate a 2D array from a 1D array, given the desired dimensions
array2D::array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2 ) 
        : arraydata( dataOneD ){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
    if (dataOneD->size() != size_dim1*size_dim2) {
        cerr << "WARNING in array2D::array2D. Inconsistent array size. ";
        cerr << "size1D=" << dataOneD->size() << ", size2D=" << size_dim1*size_dim2 
			<< "=" << size_dim1 << "*" << size_dim2 << "" << endl;
    }
}

array2D::array2D( array2D* dataTwoD )
		: arraydata(dataTwoD){
	setDim1( dataTwoD->dim1() );
	setDim2( dataTwoD->dim2() );
}

array2D::~array2D(){
}


//-----------------------------------------------------copy
void array2D::copy( const array2D& src ){
    setDim1( src.dim1() );
    setDim2( src.dim2() );
    this->arraydata::copy( src.data(), src.size() );
}

//-----------------------------------------------------get
//time-critical function, check with assert is disabled in release configuration
double array2D::get( unsigned int i, unsigned int j ) const{
	assert(i < dim1());
	assert(j < dim2());
	return arraydata::get_atIndex( j*dim1() + i );
}

//-----------------------------------------------------set
//time-critical function, check with assert is disabled in release configuration
void array2D::set( unsigned int i, unsigned int j, double value ){
	assert(i < dim1());
	assert(j < dim2());
	arraydata::set_atIndex( j*dim1() + i, value);
}


//-----------------------------------------------------getRow
int array2D::getRow( int rownum, array1D *&row ) const{
	if (rownum >= dim1() || rownum < 0){ 
		cerr << "Error in array2D::getRow. row number " << rownum << " too big or below zero." << endl; 
		return 1;
	}
	if (this->dim2()==0){
		cerr << "Error in array2D::getRow. array2D's dimension2 is zero" << endl; 
		return 2;
	}
	
	// create fresh array if 'row' doesn't have right size, otherwise, just overwrite data
	if (row->size() != this->dim2()) {
		delete row;
		try{
			row = new array1D( this->dim2() );
		}catch (std::exception& e){
			cerr << "Error in array2D::getRow. Could not allocate row." << endl;
			cerr << "Standard exception: " << e.what() << endl;
		}
	}

	//for a fixed row, i goes through columns (x-values)
	for (int i = 0; i < row->size(); i++){
		row->set( i, this->get(rownum, i) );
	}
	return 0;
}

//-----------------------------------------------------getCol
int array2D::getCol( int colnum, array1D *&col) const{
    if (colnum >= dim2() || colnum < 0){ 
		cerr << "Error in array2D::getCol. column number " << colnum << " too big or below zero." << endl; 
		return 1;
	}
	if (this->dim1()==0){
		cerr << "Error in array2D::getCol. array2D's dimension1 is zero" << endl; 
		return 2;
	}

	// create fresh array if 'col' doesn't have right size, otherwise, just overwrite data
	if (col->size() != this->dim1()) {		
    	delete col;
		try{
			col = new array1D( this->dim1() );
		}catch (std::exception& e){
			cerr << "Error in array2D::getCol. Could not allocate column." << endl; 
			cerr << "Standard exception: " << e.what() << endl;
		}
 	}
	
	//for a fixed column number, j goes through the rows (y-values)
	for (int j = 0; j < col->size(); j++){
		col->set( j, this->get(j, colnum) );
	}
	return 0;
}




	

//-----------------------------------------------------setRow
// note: there is no check for dimensions
// the array2D is updated as long as there is data in the passed array1D
// therefore, if dim(1D) < dim(2D), there will be old data, which is not overwritten
//       and, if dim(1D) > dim(2D), not all data from 1D will be present in the 2D
void array2D::setRow( int rownum, const array1D *row, int start ){
	if (!row){
		cerr << "Error in array2D::setRow. Passed 'row' not allocated." << endl;
		throw;
	}else if( start < 0 || start >= row->size() ){
		cerr << "Error in array2D::setRow. Start value " << start << " not allowed." << endl;
		throw;
	}else{
		for (int i = 0; i < row->size() && start+i < this->dim2(); i++){
			this->set( rownum, start+i, row->get(i) );
		}
	}
}

//-----------------------------------------------------setCol
void array2D::setCol( int colnum, const array1D *col, int start ){
	if (!col){
		cerr << "Error in array2D::setCol. Passed 'col' not allocated." << endl;
		throw;
	}else if( start < 0 || start >= col->size() ){
		cerr << "Error in array2D::setCol. Start value " << start << " not allowed." << endl;
		throw;
	}else{
		for (int i = 0; i < col->size() && start+i < this->dim1(); i++){
			this->set( start+i, colnum, col->get(i) );
		}
	}
}

//-----------------------------------------------------calcAvgRow
int array2D::calcAvgRow( array1D *&avgrow ) const{
	delete avgrow;
	avgrow = new array1D( this->dim1() );
	for (int c = 0; c < dim2(); c++){
		double sum = 0;
		for (int r = 0; r < dim1(); r++){
			sum += get(r,c);
		}
		avgrow->set(c, sum/dim1());
	}
	return 0;
}

//-----------------------------------------------------calcAvgCol
int array2D::calcAvgCol( array1D *&avgcol ) const{
	delete avgcol;
	avgcol = new array1D( this->dim1() );
	for (int r = 0; r < dim1(); r++){
		double sum = 0;
		for (int c = 0; c < dim2(); c++){
			sum += get(r,c);
		}
		avgcol->set(r, sum/dim2());
	}
	return 0;
}


//------------------------------------------------------------- transpose
void array2D::transpose(){
	if (this->dim1()==0 || this->dim2()==0){
		cerr << "Error in array2D::transpose. dim1=" << dim1() << ", dim2=" << dim2() << ". Nothing done." << endl;
		return;
	}
	
	array2D *old = new array2D(*this);
	
	//swap dimensions, total arraydata length stays the same
	this->setDim1(old->dim2());
	this->setDim2(old->dim1());
	for ( int i = 0; i < dim1(); i++ ){
		for ( int j = 0; j < dim2(); j++ ){
			this->set( i, j, old->get(j,i) );
		}
	}
	delete old;
}

//------------------------------------------------------------- flipud
void array2D::flipud(){
	array2D *old = new array2D(*this);
	for ( int j = 0; j < dim2(); j++ ){
		for ( int i = 0; i < dim1(); i++ ){
			this->set( i, j, old->get(dim1()-1-i,j) );
		}
	}
	delete old;
}

//------------------------------------------------------------- fliplr
void array2D::fliplr(){
	array2D *old = new array2D(*this);
	for ( int i = 0; i < dim1(); i++ ){
		for ( int j = 0; j < dim2(); j++ ){
			this->set( i, j, old->get(i,dim2()-1-i) );
		}
	}
	delete old;	
}


//------------------------------------------------------------- gradientAlongDim1
// remains constant along dim1 (columns) 
// example of gradientAlongDim1(-2,2)
//	  2  2  2  2  2
//	  1  1  1  1  1
//	  0  0  0  0  0
//	 -1 -1 -1 -1 -1
//	 -2 -2 -2 -2 -2
//
void array2D::gradientAlongDim1( double lowlim, double highlim ){	//set elements to a range of values, given by the boundaries
	for (int j = 0; j < dim2(); j++) {
		array1D *col = new array1D( dim1() );
		col->range(lowlim, highlim);
		this->setCol(j, col);
		delete col;
	}
}

//------------------------------------------------------------- gradientAlongDim2
// remains constant along dim2 (rows)
// example of gradientAlongDim2(-2,2)
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//
void array2D::gradientAlongDim2( double lowlim, double highlim ){	//set elements to a range of values, given by the boundaries
	for (int i = 0; i < dim1(); i++) {
		array1D *row = new array1D( dim2() );
		row->range(lowlim, highlim);
		this->setRow(i, row);
		delete row;
	}
}


    
//------------------------------------------------------------- getASCIIdata
std::string array2D::getASCIIdata( bool annotate ) const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "2D data has size zero." << endl;
	}else{
		if (annotate){ osst << "2D data, dim1=" << dim1() << ", dim2=" << dim2() << ", size=" << size() << endl; }
		for (int i = 0; i<dim1(); i++){
			if (annotate){ osst << " ["; }
			for (int j = 0; j<dim2(); j++){
				osst << " " << get(i, j);
			}
			if (annotate){ osst << "]"; }
			osst << endl;
		}
	}
    return osst.str();
}





//-----------------------------------------------------setters & getters
unsigned int array2D::dim1() const{
	return p_dim1;
}

void array2D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array2D::dim2() const{
	return p_dim2;
}

void array2D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}




//-----------------------------------------------------test pattern
void array2D::generateTestPattern( int type ){
    
    if (dim1()==0 || dim2()==0)
        cout << "WARNING in generateTestPattern. One or more dimensions are zero." << endl;


    cout << "Writing test pattern type " << type << ": ";
	switch (type) {
		case 0:											
			{   
                cout << "2D sinusoidal ";
				double amplitude = 20000;
				double periodX = dim1()/1;
				double periodY = dim2()/2;
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
						set(i, j, amplitude/2*(sin(2*M_PI*i/periodX)*cos(2*M_PI*j/periodY)+1) );
					}
				}
			}
            break;
		case 1:											
			{
                cout << "increment by absolute array index ";
				double val = 0;
				for (int i = 0; i < size(); i++) {
					if (val >= 65535) {
						val = 0; 
					}
					set_atIndex(i, val);
					val += 1;
				}
			}
			break;
		case 2:											
        	{
                cout << "2D centro-symmetric sine ";
				double amplitude = 20000;
				double period = dim1()/5;
                double centerX = dim1()/2;
                double centerY = dim2()/2;
                cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
                        double x = i - centerX;
                        double y = j - centerY;
                        double r = sqrt( x*x + y*y );
						set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 ) );
					}
				}
			}
            break;
		case 3:											
        	{
                cout << "2D centro-symmetric sine with circular modulation ";
				double amplitude = 20000;
				double period = dim1()/5;
                double periodMod = dim1()/10;
                double centerX = dim1()/2;
                double centerY = dim2()/2;
                cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
                        double x = i - centerX;
                        double y = j - centerY;
                        double r = sqrt( x*x + y*y );
                        double phi = atan(y/x);
						set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 + 1/10*(tan(2*M_PI*phi/periodMod)+1) ) );
					}
				}
			}
            break;
		case 4:											
        	{
                cout << "2D centro-symmetric sine with straight modulation ";
				double amplitude = 20000;
				double period = dim1()/5;
                double periodMod = dim1()/10;
                double centerX = dim1()/2;
                double centerY = dim2()/2;
                cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
                        double x = i - centerX;
                        double y = j - centerY;
                        double r = sqrt( x*x + y*y );
						set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 + cos(2*M_PI*i/periodMod)+1 ) );
					}
				}
			}
            break;
		default:                                        
            cout << "zeros";
			zeros();
			break;
	}
    cout << endl;
}




//=================================================================================
//
// CLASS IMPLEMENTATION OF array3D
//
//=================================================================================
//-----------------------------------------------------constructors & destructors
array3D::array3D( unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3 )
		: arraydata(size_dim1*size_dim2*size_dim3){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
	setDim3( size_dim3 );
}

array3D::~array3D(){
}



//-----------------------------------------------------copy
void array3D::copy( const array3D& src ){
    setDim1( src.dim1() );
    setDim2( src.dim2() );
    setDim3( src.dim3() );
    this->arraydata::copy( src.data(), src.size());
}


//-----------------------------------------------------get
//time-critical function, check with assert is disabled in release configuration
double array3D::get( unsigned int i, unsigned int j, unsigned int k ) const{
	assert(i < dim1());
	assert(j < dim2());
	assert(k < dim3());
	return arraydata::get_atIndex( k*dim1()*dim2() + j*dim1() + i );
}

//-----------------------------------------------------set
//time-critical function, check with assert is disabled in release configuration
void array3D::set( unsigned int i, unsigned int j, unsigned int k, double value ){
	assert(i < dim1());
	assert(j < dim2());
	assert(k < dim3());
	arraydata::set_atIndex( k*dim1()*dim2() + j*dim1() + i, value);
}


//-----------------------------------------------------setters & getters
unsigned int array3D::dim1() const{
	return p_dim1;
}

void array3D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array3D::dim2() const{
	return p_dim2;
}

void array3D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}

unsigned int array3D::dim3() const{
	return p_dim3;
}

void array3D::setDim3( unsigned int size_dim3 ){
	p_dim3 = size_dim3;
}


//------------------------------------------------------------- getASCIIdata
string array3D::getASCIIdata( bool annotate ) const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "3D data has size zero." << endl;
	}else{
		if (annotate){ osst << "3D data, dim1=" << dim1() << ", dim2=" << dim2() << ", dim3=" << dim3() << ", size=" << size() << endl; }
		for (int k = 0; k<dim3(); k++){
			if (annotate){ osst << " [[" << endl; }
			for (int j = 0; j<dim2(); j++){
				if (annotate){ osst << "  ["; }
				for (int i = 0; i<dim1(); i++) {
					osst << " " << get(i, j, k);
				}//k
				if (annotate){ osst << "]"; }
				osst << endl;
			}//j
			if (annotate){ osst << "]]"; }
			osst << endl;
		}//i
	}
    return osst.str();
}




//=================================================================================
//
// CLASS IMPLEMENTATION OF array4D
//
//=================================================================================
//-----------------------------------------------------constructors & destructors
array4D::array4D( unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3, unsigned int size_dim4 )
		: arraydata(size_dim1*size_dim2*size_dim3*size_dim4){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
	setDim3( size_dim3 );
	setDim4( size_dim4 );
}

array4D::array4D( array1D *dataOneD, 
		unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3, unsigned int size_dim4 )
		: arraydata( dataOneD ){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
	setDim3( size_dim3 );
	setDim4( size_dim4 );
	if (dataOneD->size() != size_dim1*size_dim2*size_dim3*size_dim4) {
        cerr << "WARNING in array4D::array4D. Inconsistent array size. ";
        cerr << "size1D=" << dataOneD->size() 
			<< ", size4D=" << size_dim1*size_dim2*size_dim3*size_dim4 << endl;
    }
}


array4D::~array4D(){
}



//-----------------------------------------------------copy
void array4D::copy( const array4D& src ){
    setDim1( src.dim1() );
    setDim2( src.dim2() );
    setDim3( src.dim3() );
	setDim4( src.dim4() );
    this->arraydata::copy( src.data(), src.size() );
}


//-----------------------------------------------------get
//time-critical function, check with assert is disabled in release configuration
double array4D::get( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const{
	assert(i < dim1());
	assert(j < dim2());
	assert(k < dim3());
	assert(l < dim4());
	return arraydata::get_atIndex( l*dim1()*dim2()*dim3() + k*dim1()*dim2() + j*dim1() + i );
}

//-----------------------------------------------------set
//time-critical function, check with assert is disabled in release configuration
void array4D::set( unsigned int i, unsigned int j, unsigned int k, unsigned int l, double value ){
	assert(i < dim1());
	assert(j < dim2());
	assert(k < dim3());
	assert(l < dim4());
	arraydata::set_atIndex( l*dim1()*dim2()*dim3() + k*dim1()*dim2() + j*dim1() + i, value);
}


//-----------------------------------------------------setters & getters
unsigned int array4D::dim1() const{
	return p_dim1;
}

void array4D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array4D::dim2() const{
	return p_dim2;
}

void array4D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}

unsigned int array4D::dim3() const{
	return p_dim3;
}

void array4D::setDim3( unsigned int size_dim3 ){
	p_dim3 = size_dim3;
}

unsigned int array4D::dim4() const{
	return p_dim4;
}

void array4D::setDim4( unsigned int size_dim4 ){
	p_dim4 = size_dim4;
}


//------------------------------------------------------------- getASCIIdata
string array4D::getASCIIdata( bool annotate ) const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "4D data has size zero." << endl;
	}else{
		if (annotate){ osst << "4D data, dim1=" << dim1() << ", dim2=" << dim2() << ", dim3=" << dim3() << ", dim4=" << dim4() << ", size=" << size() << endl; }
		for (int l = 0; l<dim4(); l++){
			for (int k = 0; k<dim3(); k++){
				if (annotate){ osst << " [[" << endl; }
				for (int j = 0; j<dim2(); j++){
					if (annotate){ osst << "  ["; }
					for (int i = 0; i<dim1(); i++) {
						osst << " " << get(i, j, k, l);
					}//i
					if (annotate){ osst << "]"; }
					osst << endl;
				}//j
				if (annotate){ osst << "]]"; }
				osst << endl;
			}//k
		}//l
	}
    return osst.str();
}




//------------------------------------------------------------- 
// mainly for creating 'raw' CSPAD images right now, where
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
//
void array4D::getRepresentationIn2D( array2D *&img ){
	delete img;
	img = new array2D(dim1()*dim4(), dim2()*dim3());	
//	cout << "4D data, dim1=" << dim1() << ", dim2=" << dim2() << ", dim3=" << dim3() << ", dim4=" << dim4() << ", size=" << size() << endl;
//	cout << "2D data, dim1=" << img->dim1() << ", dim2=" << img->dim2() << endl; 
	for (int q = 0; q<dim4(); q++){
		for (int s = 0; s<dim3(); s++){
			for (int c = 0; c<dim2(); c++){
				for (int r = 0; r<dim1(); r++) {
					int row = q*dim1() + r;
					int col = s*dim2() + c;
					img->set( row, col, get(r, c, s, q) );
				}//i
			}//j
		}//k
	}//l	
	
	//transpose to be conform with cheetah's convention
	img->transpose();
}

