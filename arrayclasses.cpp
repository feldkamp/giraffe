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

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <fstream>
using std::ofstream;

#include <cmath>
#include <cassert>




//=================================================================================
//
// CLASS IMPLEMENTATION OF arraydata
//
//=================================================================================
//arraydata::arraydata(){
//    init();
//	p_data = new double[0];
//}

arraydata::arraydata( unsigned int size_val ){
    init();
	p_size = size_val;
	p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;

	zeros();					// set all elements to zero initially
}

arraydata::arraydata( const int16_t *CArray, const unsigned int size_val ){
    init();
    p_size = size_val;
    p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;
    
    //fill array with data, convert int to double before
    for (int i = 0; i < size_val; i++) {
        p_data[i] = (double)CArray[i];
    }
}

arraydata::arraydata( const float *CArray, const unsigned int size_val ) {
    init();
    p_size = size_val;
    p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;
    
    //fill array with data, convert float to double before
    for (int i = 0; i < size_val; i++) {
        p_data[i] = (double)CArray[i];
    }
}

arraydata::arraydata( const arraydata &src ){                       //copy constructor
    init();
    copy( src );
}

arraydata & arraydata::operator=(const arraydata & src){
    if ( this != &src ){
		cout << "DEBUG: ASSIGNMENT OPERATOR FOR arraydata" << endl;
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

//-----------------------------------------------------copy
void arraydata::copy( const double* src_data, const unsigned int arraysize ){		//src type: double c-array
    p_size = arraysize;
    if (p_size > 0){
		this->destroy();
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = src_data[i];
        }
    }else{
        p_data = NULL;
    }
}

void arraydata::copy( const float* src_data, const unsigned int arraysize ){		//src type: double c-array
    p_size = arraysize;
    if (p_size > 0){
		this->destroy();
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = (double)src_data[i];
        }
    }else{
        p_data = NULL;
    }
}

void arraydata::copy( const int* src_data, const unsigned int arraysize ){		//src type: int c-array
    p_size = arraysize;
    if (p_size > 0){
		this->destroy();
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = (double) src_data[i];	//convert to double
        }
    }else{
        p_data = NULL;
    }
}

void arraydata::copy( const arraydata& src ){						//src type: arraydata object
    this->copy( src.data(), src.size());
}


double *arraydata::data() const{
    return (double *)p_data;
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

//--------------------------------------------------------------------
double arraydata::calcMin() const{
	double tempmin = +INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atIndex(i) < tempmin) {
			tempmin = get_atIndex(i);
		}
	}
	return tempmin;
}

//--------------------------------------------------------------------
double arraydata::calcMax() const{
	double tempmax = -INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atIndex(i) > tempmax) {
			tempmax = get_atIndex(i);
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

//------------------------------------------------------------- getASCIIdata
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

//------------------------------------------------------------- getASCIIdata
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

array1D::array1D( int16_t *CArray, unsigned int size_dim1 )
        : arraydata( CArray, size_dim1 ){
    setDim1( size_dim1 );
}

array1D::array1D( float *CArray, unsigned int size_dim1 )
		: arraydata( CArray, size_dim1 ) {
    setDim1( size_dim1 );
}


//constructor to generate a 1D array from a 2D array
array1D::array1D( array2D* dataTwoD ) 
        : arraydata( dataTwoD->size() ){
 	setDim1( dataTwoD->size() );

    //copy contents of dataTwoD to this array1D object
    p_size = dataTwoD->size();
    if (p_size > 0){
        for (int i = 0; i < p_size; i++) {
            p_data[i] = dataTwoD->get_atIndex(i);
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
string array1D::getASCIIdata() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "1D data has size zero." << endl;
	}else{
		osst << "1D data, dim1=" << dim1() << ", size=" << size() << endl;
		osst << " [";
		for (int i = 0; i<dim1(); i++) {
			osst << " " << get(i);
		}
		osst << "]" << endl;
	}
    return osst.str();
}


//------------------------------------------------------------- writeToASCII
int array1D::writeToASCII( std::string filename ) const{
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
	return 0;
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


//constructor to generate a 2D array from a 1D array, given the desired dimensions
array2D::array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2) 
        : arraydata( size_dim1*size_dim2 ){
    
    if (!dataOneD) {
        cerr << "WARNING in array2D::array2D. Input array1D was not allocated! Nothing copied!" << endl;
        return;
    }
    if (dataOneD->size() != size_dim1*size_dim2) {
        cerr << "WARNING in array2D::array2D. Inconsistent array size. ";
        cerr << "size1D=" << dataOneD->size() << ", size2D=" << size_dim1*size_dim2 
			<< "=" << size_dim1 << "*" << size_dim2 << "" << endl;
    }
    
 	setDim1( size_dim1 );
	setDim2( size_dim2 );

    //copy contents of dataOneD to this array2D
    p_size = dataOneD->size();
    if (p_size > 0){
        for (int i = 0; i < p_size; i++) {
            p_data[i] = dataOneD->get_atIndex(i);
        }
    }else{
        p_data = NULL;
    }
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


//-----------------------------------------------------more data accessors
int array2D::getRow( int rownum, array1D *&row ) const{
	if (rownum >= dim1() || rownum < 0){ 
		cerr << "Error in array2D::getRow. row number " << rownum << " too big or below zero." << endl; 
		return 1;
	}
	if (this->dim2()==0){
		cerr << "Error in array2D::getRow. array2D's dimension2 is zero" << endl; 
		return 2;
	}
	
	// create new array if 'row' doesn't have right size, otherwise, just overwrite data
	if (row->size() != this->dim2()) {
		delete row;
		row = new array1D( this->dim2() );
		if (!row){ 
			cerr << "Error in array2D::getRow. Could not allocate row." << endl;
			return 3;
		}
	}

	//for a fixed row, i goes through columns (x-values)
	for (int i = 0; i < row->size(); i++){
		row->set( i, this->get(rownum, i) );
	}
	return 0;
}


int array2D::getCol( int colnum, array1D *&col) const{
    if (colnum >= dim2() || colnum < 0){ 
		cerr << "Error in array2D::getCol. column number " << colnum << " too big or below zero." << endl; 
		return 1;
	}
	if (this->dim1()==0){
		cerr << "Error in array2D::getCol. array2D's dimension1 is zero" << endl; 
		return 2;
	}

	// create new array if 'col' doesn't have right size, otherwise, just overwrite data
	if (col->size() != this->dim1()) {		
    	delete col;
    	col = new array1D( this->dim1() );
		if (!col){ 
			cerr << "Error in array2D::getCol. Could not allocate column." << endl; 
			return 3;
		}
 	}
	
	//for a fixed column number, j goes through the rows (y-values)
	for (int j = 0; j < col->size(); j++){
		col->set( j, this->get(j, colnum) );
	}
	return 0;
}


//-----------------------------------------------------setRow/setCol
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

//------------------------------------------------------------- transpose
void array2D::transpose(){
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
std::string array2D::getASCIIdata() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "2D data has size zero." << endl;
	}else{
		osst << "2D data, dim1=" << dim1() << ", dim2=" << dim2() << ", size=" << size() << endl;
		for (int i = 0; i<dim1(); i++){
			osst << " [";
			for (int j = 0; j<dim2(); j++){
				osst << " " << get(i, j);
			}
			osst << "]" << endl;
		}
	}
    return osst.str();
}



//------------------------------------------------------------- writeToASCII
int array2D::writeToASCII( std::string filename, int format ) const{
	ofstream fout( filename.c_str() );
	switch (format){
		case 1:
			fout << arraydata::getASCIIdataAsColumn();
		case 2:
			fout << arraydata::getASCIIdataAsRow();
		default:
			fout << getASCIIdata();
	}
	fout.close();
	return 0;
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
string array3D::getASCIIdata() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "3D data has size zero." << endl;
	}else{
		osst << "3D data, dim1=" << dim1() << ", dim2=" << dim2() << ", dim3=" << dim3() << ", size=" << size() << endl;
		for (int k = 0; k<dim3(); k++){
			osst << " [[" << endl;
			for (int j = 0; j<dim2(); j++){
				osst << "  [";
				for (int i = 0; i<dim1(); i++) {
					osst << " " << get(i, j, k);
				}//k
				osst << "]" << endl;
			}//j
			osst << "]]" << endl;
		}//i
	}
    return osst.str();
}



//------------------------------------------------------------- writeToASCII
int array3D::writeToASCII( std::string filename ) const{
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
    return 0;
}




