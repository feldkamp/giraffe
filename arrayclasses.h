/*
 *  ArrayClasses.h
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 04/27/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */

#ifndef _arrayclasses_h
#define _arrayclasses_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <cassert>
#include <exception>

//=================================================================================
//
// arraydata (the base class for array1D, array2D, array3D)
//
//=================================================================================
template<typename T>
class arraydata {

private:
	std::vector<T> p_data;
	
public:
	//-----------------------------------------------------------------------------------------constructors
	arraydata( const unsigned int sizeval = 1 ){
		p_data.assign(sizeval, 0);						//assign the value zero to a vector of length sizeval
	}
	arraydata( const arraydata *array ){
		copy( array );
	}
	arraydata( const std::vector<T> vec ){
		copy( vec );
	}
	template <typename arrayT> arraydata( const arrayT *CArray, const unsigned int size_val ){
		copy( CArray, size_val );
	}
	arraydata( const arraydata &src ){                       	//copy constructor
    	copy( src );
	} 
			
	//-----------------------------------------------------------------------------------------copy functions
	void copy( const arraydata *src ){
		p_data = src->data();
	}
	void copy( const arraydata& src ){
		p_data = src.data();
	}
	void copy( const std::vector<T> vec ){
		p_data = vec;
	}

	// templated function for c-style arrays of basic types (float, double, int, ...)
	template <typename arrayT> void copy( const arrayT *src, const unsigned int arraysize ){
		p_data.clear();
		p_data.reserve( arraysize );
		if (arraysize > 0 && src){
			
			//copy all elements from src to internal data object
			//casting them to type T
			for (unsigned int i = 0; i < arraysize; i++) {
				p_data.push_back( (T)(src[i]) );
			}
		}else{
			std::cerr << "Error in arraydata::copy! src=" << src << ", size=" << arraysize << std::endl;
		}
	}
	
								 
	//-----------------------------------------------------------------------------------------operators
    arraydata & operator=(const arraydata & src){               //assignment operator
		if ( this != &src ){
			copy( src );
		}
		return *this;
	}
	
	T const& operator[]( int i ) const	{ return p_data[i]; }
	T& operator[]( int i )				{ return p_data[i]; }	

	std::vector<T> data() const 		{ return p_data; }		//return internal raw data
    

	//-----------------------------------------------------------------------------------------data accessors
	// 'atIndex' functions:
	// the following functions change or return properties
	// assuming a one-dimensional array (consistent with the internal data structure), 
    // no matter what dimension the actual subclass may have
	unsigned int size() const							{ return (int)p_data.size(); }	//total size of the array	
	T get_atIndex( unsigned int i) const				{ return p_data[i]; }		// get element value
	void set_atIndex( unsigned int i, T val)			{ p_data[i] = val; }		// set element value
    
	//-----------------------------------------------------------------------------------------basic assignments
    void zeros(){												//set all elements to 0
		for (unsigned int i = 0; i < size(); i++) {
			set_atIndex(i, 0);
		}
	}
	
    void zeros( unsigned int start, unsigned int stop ){		//set all elements between start and stop to 0 (incl. start, excl. stop)
		if ( start >= stop || stop > size() ){
			std::cerr << "Error in arraydata::zero("<< start << ", " << stop << "). Check boundaries." << std::endl;
			throw;
		}
		for (int i = start; i < stop; i++) {
			set_atIndex(i, 0);
		}
	}	
	
	void ones(){												//set all elements to 1
		for (unsigned int i = 0; i < size(); i++) {
			set_atIndex(i, 1);
		}	
	}
	
    void ones( unsigned int start, unsigned int stop ){			//set all elements between start and stop to 1 (incl. start, excl. stop)
	if ( start >= stop || stop > size() ){
			std::cerr << "Error in arraydata::ones("<< start << ", " << stop << "). Check boundaries." << std::endl;
			throw;
		}
		for (int i = start; i < stop; i++) {
			set_atIndex(i, 1);
		}
	}
	
	void range( double neg, double pos ){						//set elements to a range of values, given by the boundaries
		double delta = (pos-neg)/(size()-1);
		for (unsigned int i = 0; i < size(); i++) {
			set_atIndex(i, neg+delta*i);
		}	
	}
	//-------------------------------------------------------------------- calcMin
	// ignores the position of the minimum
	T calcMin() const{
		int pos = 0;					
		return this->calcMin(pos);
	}

	//-------------------------------------------------------------------- calcMin
	// returns the position of the minimum
	T calcMin(int &pos) const{
		double tempmin = +INFINITY;
		for (unsigned int i = 0; i < size(); i++) {
			if (get_atIndex(i) < tempmin) {
				if ( get_atIndex(i) != -INFINITY ){
					tempmin = get_atIndex(i);
					pos = i;
				}else{
					std::cout << "WARNING in arraydata::calcMax(). Data contains '-infinity' at index " << i << ". Skipping this value." << std::endl;
				}
			}
		}
		return tempmin;
	}

	//-------------------------------------------------------------------- calcMax
	// ignores the position of the minimum
	T calcMax() const{
		int pos = 0;					
		return this->calcMax(pos);
	}

	//-------------------------------------------------------------------- calcMax
	// returns the position of the maximum
	T calcMax(int &pos) const{
		double tempmax = -INFINITY;
		for (unsigned int i = 0; i < size(); i++) {
			if (get_atIndex(i) > tempmax) {
				if ( get_atIndex(i) != INFINITY ){
					tempmax = get_atIndex(i);
					pos = i;
				}else{
					std::cout << "WARNING in arraydata::calcMax(). Data contains 'infinity' at index " << i << ". Skipping this value." << std::endl;
				}
			}
		}
		return tempmax;
	}

	//------------------------------------------------------------- calcSum
	T calcSum() const{
		double sum = 0.;
		for (unsigned int i = 0; i < size(); i++) {
			sum += get_atIndex(i);
		}	
		return sum;
	}

	//------------------------------------------------------------- calcAvg
	double calcAvg() const{
		double avg = calcSum() / ((double)size());
		return avg;
	}

	//------------------------------------------------------------- calcStDev
	double calcStDev() const{
		double avg = 0.;	//value is effectively discarded
		return calcStDev(avg);
	}

	double calcStDev( double &avg ) const{
		double sum_sqrd = 0.;
		avg = this->calcAvg();
		for (unsigned int i = 0; i < size(); i++){
			double x = this->get_atIndex(i) - avg;
			sum_sqrd += x*x;
		}
		double sigma = sqrt( sum_sqrd/((double)size()) );
		return sigma;
	}

	//------------------------------------------------------------- getASCIIdataAsRow
	std::string getASCIIdataAsRow() const{
		std::ostringstream osst;
		if (size() == 0) {
			osst << "Array data has size zero." << std::endl;
		}else{
			for (unsigned int i = 0; i<size(); i++) {
				osst << get_atIndex(i) << " ";
			}
			osst << std::endl;
		}
		return osst.str();
	}

	//------------------------------------------------------------- getASCIIdataAsColumn
	std::string getASCIIdataAsColumn() const{
		std::ostringstream osst;
		if (size() == 0) {
			osst << "Array data has size zero." << std::endl;
		}else{
			for (unsigned int i = 0; i<size(); i++) {
				osst << get_atIndex(i) << std::endl;
			}
		}
		return osst.str();
	}


	//--------------------------------------------------------------------array math
	//multiply each element by a numerical factor
	int addValue( double val ){
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i) + val);
		}
		return 0;
	}

	//multiply each element by a numerical factor
	int subtractValue( double val ){
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i) - val);
		}
		return 0;
	}

	//multiply each element by a numerical value
	int multiplyByValue( double value ){
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i) * value);
		}
		return 0;
	}

	//divide each element by a numerical value (enforcing float division)
	int divideByValue( double value ){
		if (value != 0){
			for (unsigned int i = 0; i<this->size(); i++) {
				this->set_atIndex(i, ((double) this->get_atIndex(i)) / (double)value);
			}
		}else{
			std::cerr << "Error in arraydata::divideByValue. Division by zero requested! Nothing done." << std::endl;
			return 1;
		}
		return 0;
	}

	int addValue_atIndex( unsigned int i, double value ){
		set_atIndex( i, get_atIndex(i) + value );
		return 0;
	}

	int subtractValue_atIndex( unsigned int i, double value ){
		set_atIndex( i, get_atIndex(i) - value );
		return 0;
	}

	int multiplyByValue_atIndex( unsigned int i, double value ){
		set_atIndex( i, get_atIndex(i) * value );
		return 0;
	}

	int divideByValue_atIndex( unsigned int i, double value ){
		if (value != 0){
			set_atIndex( i, get_atIndex(i) / (double)value );
		}else{
			std::cerr << "Error in arraydata::divideByValue_atIndex("<< i << "). Division by zero requested! Nothing done." << std::endl;
			return 1;		
		}
		return 0;
	}

	int increment_atIndex( unsigned int i ){
		set_atIndex( i, get_atIndex(i) - 1 );
		return 0;
	}

	int decrement_atIndex( unsigned int i ){
		set_atIndex( i, get_atIndex(i) + 1 );
		return 0;
	}



	//add each element by an element from a second array
	int addArrayElementwise( const arraydata *secondArray ){
		if (!secondArray){
			std::cerr << "Error in arraydata::addArrayElementwise! Second array not allocated. ";
			return 2;
		}
		if (this->size() != secondArray->size()){
			std::cerr << "Error in arraydata::addArrayElementwise! Array sizes don't match. ";
			std::cerr << "(this array size " << this->size() << " != second array size " << secondArray->size() << "). Operation not performed."<< std::endl;
			return 1;
		}
		
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i)+secondArray->get_atIndex(i));
		}
		return 0;
	}

	//subtract each element by an element from a second array
	int subtractArrayElementwise( const arraydata *secondArray ){
		if (!secondArray){
			std::cerr << "Error in arraydata::subtractArrayElementwise! Second array not allocated. ";
			return 2;
		}
		if (this->size() != secondArray->size()){
			std::cerr << "Error in arraydata::subtractArrayElementwise! Array sizes don't match. ";
			std::cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< std::endl;
			return 1;
		}
		
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i)-secondArray->get_atIndex(i));
		}
		return 0;
	}

	//multiply each element by an element from a second array
	int multiplyByArrayElementwise( const arraydata *secondArray ){
		if (!secondArray){
			std::cerr << "Error in arraydata::multiplyArrayElementwise! Second array not allocated. ";
			return 2;
		}
		if (this->size() != secondArray->size()){
			std::cerr << "Error in arraydata::multiplyArrayElementwise! Array sizes don't match. ";
			std::cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< std::endl;
			return 1;
		}
		
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i) * secondArray->get_atIndex(i) );
		}
		return 0;
	}

	//divide each element by an element from a second array
	int divideByArrayElementwise( const arraydata *secondArray ){
		if (!secondArray){
			std::cerr << "Error in arraydata::divideArrayElementwise! Second array not allocated. ";
			return 2;
		}
		if (this->size() != secondArray->size()){
			std::cerr << "Error in arraydata::divideArrayElementwise! Array sizes don't match. ";
			std::cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< std::endl;
			return 1;
		}
		
		for (unsigned int i = 0; i<this->size(); i++) {
			this->set_atIndex(i, this->get_atIndex(i)/(double)secondArray->get_atIndex(i));
		}
		return 0;
	}


	//------------------------------------------------------------- applyMask
	int applyMask( arraydata* mask, double checkval = 0.5, double rejectval = 0. ){
		if (size() != mask->size()){
			std::cerr << "Error in arraydata::applyMask! Data size " << size() << " doesn't match mask size " << mask->size() << "." << std::endl;
			return 1;
		}
		for (unsigned int i = 0; i < size(); i++){
			if (mask->get_atIndex(i) < checkval){
				//replace by fill value
				set_atIndex(i, rejectval);
			}else{
				//do nothing, keep the original data
			}
		}
		return 0;
	}

	//------------------------------------------------------------- histogram
	int getHistogramInBoundaries( std::vector<int> &hist, std::vector<double> &bins, unsigned int nBins, double min, double max ){
		double range = max-min;
		double binwidth = range/(double)(nBins);
		double bindelta = range/(double)(nBins-1);
		//cout << "min: " << min << ", max:" << max << ", range: " << range << ", binwidth: " << binwidth << ", nBins: " << nBins << endl;
		
		if (range == 0 ){
			std::cerr << "WARNING in arraydata::getHistogramInBoundaries. max == min == " << max << std::endl;
		}
		
		//make bins array
		bins.assign(nBins, 0);
		for (unsigned int b = 0; b < nBins; b++){
			bins[b] = min + b*binwidth;
		}
		
		//make histogram
		hist.assign(nBins, 0);
		for (unsigned int i = 0; i < this->size(); i++){
			unsigned int binnum = (unsigned int) floor( (this->get_atIndex(i)-min) / (double)bindelta );
			if (binnum < nBins){
				hist[binnum]++;
			}
		}

		//find maximum of histogram and return that value		
		std::vector<int>::iterator it_max = std::max_element( hist.begin(), hist.end() );
		int histmax = *it_max;
		return histmax;
	}


	std::string getHistogramInBoundariesASCII( unsigned int nBins, double min, double max ){
		std::vector<int> hist;
		std::vector<double> bins;
		getHistogramInBoundaries( hist, bins, nBins, min, max );

		std::vector<int>::iterator it_max = std::max_element( hist.begin(), hist.end() );
		int histmax = *it_max;
		int maxpos = (int) std::distance( hist.begin(), it_max );
		int nMarkers = 20;
		int numberPerMarker = (int) floor(histmax/nMarkers);
		if (numberPerMarker < 1){
			numberPerMarker = 1;
		}
		
		double binsize = (max-min)/(double)(nBins+1);
		std::ostringstream osst;
		osst << std::setw(8) << "#" << " (" << std::setw(12) << "low" << ", " 
			<< std::setw(12) << "high" << "): " << std::setw(12) << "occurrence" << std::endl;
		osst << "---------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i < bins.size(); i++){
			osst << std::setw(8) << i << " (" << std::setw(12) << bins[i] << ", " 
				<< std::setw(12) << bins[i]+binsize << "): " 
				<< std::setw(12) << hist[i] << "  |";
			for (int m = 0; m < (hist[i]/(double)numberPerMarker); m++){ osst << "*"; }
			osst << std::endl;
		}
		osst << "---------------------------------------------------" << std::endl;
		osst << "max:" << std::setw(4) << maxpos << " (" << std::setw(12) << bins[maxpos] << ", " 
			<< std::setw(12) << bins[maxpos]+binsize << "): " << std::setw(12) << hist[maxpos] << "  |" << std::endl;
				
		return osst.str();
	}

	int getHistogram( std::vector<double> &hist, std::vector<double> &bins, unsigned int nBins = 20 ){
		double min = this->calcMin();
		double max = this->calcMax();
		return this->getHistogramInBoundaries( hist, bins, nBins, min, max );
	}

	std::string getHistogramASCII( unsigned int nBins = 20 ){
		double min = this->calcMin();
		double max = this->calcMax();
		return this->getHistogramInBoundariesASCII( nBins, min, max );
	}

};




//=================================================================================
//
// array1D -- 1-dimensional array
//
//=================================================================================
template<typename T>
class array1D : public arraydata<T>{
	
private:
	int p_dim1;
	
	void setDim( unsigned int dim1 ){
		p_dim1 = dim1;
	}
			
public:
	//-----------------------------------------------------------------------------------------constructors
	array1D( unsigned int size_dim1 = 1) 
			: arraydata<T>(size_dim1){
		setDim( size_dim1 );
	}                          
    array1D( arraydata<T>* array )
			: arraydata<T>( array ) {
		setDim( array->size() );
	}
	array1D( std::vector<T> vec )
			: arraydata<T>( vec ) {
		setDim( (int) vec.size() );
	}
	template <typename arrayT> array1D( arrayT *CArray, unsigned int size_dim1 )
        	: arraydata<T>( CArray, size_dim1 ){
		setDim( size_dim1 );
	}

	//-----------------------------------------------------------------------------------------copy functions    
    void copy( const array1D& src ){
		setDim( src.size() );
		arraydata<T>::copy( src );
	}
	
	unsigned int dim1() const 						{ return p_dim1; }
	T get( unsigned int i ) const					{ return arraydata<T>::get_atIndex(i); }
	void set( unsigned int i, T value )				{ arraydata<T>::set_atIndex(i, value); }
    
    std::string getASCIIdata( bool annotate=1 ) const{
		std::ostringstream osst;
		if (this->size() == 0) {
			osst << "1D data has size zero." << std::endl;
		}else{
			if (annotate){
				osst << "1D data, dim1=" << dim1() << ", size=" << this->size() << std::endl;
				osst << " [";
			}
			for (int i = 0; i<dim1(); i++) {
				osst << " " << get(i);
			}
			if (annotate){osst << "]";}
			osst << std::endl;
		}
		return osst.str();
	}
};




//=================================================================================
//
// array2D -- 2-dimensional array
//
// this class, while as generic as possible, adopts the convention (rows, columns)
// of packages like matlab or python for its arguments
//=================================================================================
template<typename T>
class array2D : public arraydata<T>{
private:
	int p_dim1;
	int p_dim2;

	void setDim( unsigned int dim1, unsigned int dim2 ){
		p_dim1 = dim1;
		p_dim2 = dim2;
	}
			
public:
	//-----------------------------------------------------------------------------------------constructors
	//empty 2D array 
	array2D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1 )
			: arraydata<T>(size_dim1*size_dim2){
		setDim( size_dim1, size_dim2 );
	}
	
	//init with copy of input 2D array
	array2D( array2D<T>* array ) 
			: arraydata<T>(array){
		setDim( array->dim1(), array->dim2() );
	}
	
	//init with any type of pointer to array (and two given dimensions), intended for classic C-arrays
	template <typename pointerT> 
	array2D( pointerT *dataCArray, unsigned int size_dim1, unsigned int size_dim2 )
			: arraydata<T>( dataCArray, size_dim1*size_dim2 ){
		setDim( size_dim1, size_dim2 );
	}

	//-----------------------------------------------------------------------------------------copy functions
    void copy( const array2D<T>& src ){
		setDim( src.dim1(), src.dim2() );
		arraydata<T>::copy( src.data() );
	}
	
	//-----------------------------------------------------------------------------------------data accessors
	unsigned int dim1() const		{ return p_dim1; }
	unsigned int dim2() const		{ return p_dim2; }
	
	T get( unsigned int i, unsigned int j ) const{
		assert(i < dim1());
		assert(j < dim2());
		return arraydata<T>::get_atIndex( j*dim1() + i );
	}
	
	void set( unsigned int i, unsigned int j, T value ){
		assert(i < dim1());
		assert(j < dim2());
		arraydata<T>::set_atIndex( j*dim1() + i, value);
	}

	//returns the arraydata index corresponding to the pair i, j
	unsigned int arrayIndex( unsigned int i, unsigned int j ) const {
		return (unsigned int)(j * dim1() + i);
	}

	
	//-----------------------------------------------------------------------------------------data extraction
	//returns one dimensional column, extracted at the specified column number
    int getCol( int colnum, array1D<double> *&col ) const						
	{
		if (colnum >= dim2() || colnum < 0){ 
			std::cerr << "Error in array2D::getCol. column number " << colnum << " too big or below zero." << std::endl; 
			return 1;
		}
		if (this->dim1()==0){
			std::cerr << "Error in array2D::getCol. array2D's dimension1 is zero" << std::endl; 
			return 2;
		}

		// create fresh array if 'col' doesn't have right size, otherwise, just overwrite data
		if (col->size() != this->dim1()) {		
			delete col;
			try{
				col = new array1D<double>( dim1() );
			}catch (std::exception& e){
				std::cerr << "Error in array2D::getCol. Could not allocate column." << std::endl; 
				std::cerr << "Standard exception: " << e.what() << std::endl;
			}
		}
		
		//for a fixed column number, j goes through the rows (y-values)
		for (unsigned int j = 0; j < col->size(); j++){
			col->set( j, this->get(j, colnum) );
		}
		return 0;
	}
	
	//sets a one-dimensional column, beginning at a 'start' value
	void setCol( int colnum, const array1D<double> *col, int start=0 )				
	{
		if (!col){
			std::cerr << "Error in array2D::setCol. Passed 'col' not allocated." << std::endl;
			throw;
		}else if( start < 0 || start >= col->size() ){
			std::cerr << "Error in array2D::setCol. Start value " << start << " not allowed." << std::endl;
			throw;
		}else{
			for (unsigned int i = 0; i < col->size() && start+i < this->dim1(); i++){
				this->set( start+i, colnum, col->get(i) );
			}
		}
	}		
    
	//returns one-dimensional 'row' or 'col'
    int getRow( int rownum, array1D<double> *&row ) const 	                   
	{
		if (rownum >= dim1() || rownum < 0){ 
			std::cerr << "Error in array2D::getRow. row number " << rownum << " too big or below zero." << std::endl; 
			return 1;
		}
		if (this->dim2()==0){
			std::cerr << "Error in array2D::getRow. array2D's dimension2 is zero" << std::endl; 
			return 2;
		}
		
		// create fresh array if 'row' doesn't have right size, otherwise, just overwrite data
		if (row->size() != this->dim2()) {
			delete row;
			try{
				row = new array1D<double>( this->dim2() );
			}catch (std::exception& e){
				std::cerr << "Error in array2D::getRow. Could not allocate row." << std::endl;
				std::cerr << "Standard exception: " << e.what() << std::endl;
			}
		}

		//for a fixed row, i goes through columns (x-values)
		for (unsigned int i = 0; i < row->size(); i++){
			row->set( i, this->get(rownum, i) );
		}
		return 0;
	}
	
	void setRow( int rownum, const array1D<double> *row, int start=0 ){
		if (!row){
			std::cerr << "Error in array2D::setRow. Passed 'row' not allocated." << std::endl;
			throw;
		}else if( start < 0 || start >= row->size() ){
			std::cerr << "Error in array2D::setRow. Start value " << start << " not allowed." << std::endl;
			throw;
		}else{
			for (int i = 0; i < row->size() && start+i < this->dim2(); i++){
				this->set( rownum, start+i, row->get(i) );
			}
		}
	}
	
	int calcAvgRow( array1D<double> *&row ) const{
		delete row;
		row = new array1D<double>( this->dim1() );
		for (int c = 0; c < dim2(); c++){
			double sum = 0;
			for (int r = 0; r < dim1(); r++){
				sum += get(r,c);
			}
			row->set(c, sum/dim1());
		}
		return 0;
	}
	
	int calcAvgCol( array1D<double> *&col ) const{
		delete col;
		col = new array1D<double>( this->dim1() );
		for (unsigned int r = 0; r < dim1(); r++){
			double sum = 0;
			for (unsigned int c = 0; c < dim2(); c++){
				sum += get(r,c);
			}
			col->set(r, sum/dim2());
		}
		return 0;
	}
	
	void transpose()													//transpose (dim1,dim2) --> (dim2,dim1)
	{
		if (this->dim1()==0 || this->dim2()==0){
			std::cerr << "Error in array2D::transpose. dim1=" << dim1() << ", dim2=" << dim2() << ". Nothing done." << std::endl;
			return;
		}
		
		array2D<double> *old = new array2D<double>(*this);
		
		//swap dimensions, total arraydata length stays the same
		setDim( old->dim2(),old->dim1() );
		for ( int i = 0; i < dim1(); i++ ){
			for ( int j = 0; j < dim2(); j++ ){
				this->set( i, j, old->get(j,i) );
			}
		}
		delete old;
	}
	
	void flipud()														//flip up-down
	{
		array2D<double> *old = new array2D<double>(*this);
		for ( int j = 0; j < dim2(); j++ ){
			for ( int i = 0; i < dim1(); i++ ){
				this->set( i, j, old->get(dim1()-1-i,j) );
			}
		}
		delete old;
	}
	
	void fliplr()														//flip left-right
	{
		array2D<double> *old = new array2D<double>(*this);
		for ( int i = 0; i < dim1(); i++ ){
			for ( int j = 0; j < dim2(); j++ ){
				this->set( i, j, old->get(i,dim2()-1-i) );
			}
		}
		delete old;	
	}
	
	
	//linear gradient along one dimension, same along other dimension
	void gradientAlongDim1( double lowlim, double highlim ){	//set elements to a range of values, given by the boundaries
		for (int j = 0; j < dim2(); j++) {
			array1D<double> *col = new array1D<double>( dim1() );
			col->range(lowlim, highlim);
			this->setCol(j, col);
			delete col;
		}
	}
				
	void gradientAlongDim2( double lowlim, double highlim ){	//set elements to a range of values, given by the boundaries
		for (int i = 0; i < dim1(); i++) {
			array1D<double> *row = new array1D<double>( dim2() );
			row->range(lowlim, highlim);
			this->setRow(i, row);
			delete row;
		}
	}
		
	
    std::string getASCIIdata( bool annotate=1 ) const{
		std::ostringstream osst;
		if (this->size() == 0) {
			osst << "2D data has size zero." << std::endl;
		}else{
			if (annotate){ osst << "2D data, dim1=" << dim1() << ", dim2=" << dim2() << ", size=" << this->size() << std::endl; }
			for (int i = 0; i<dim1(); i++){
				if (annotate){ osst << " ["; }
				for (int j = 0; j<dim2(); j++){
					osst << " " << get(i, j);
				}
				if (annotate){ osst << "]"; }
				osst << std::endl;
			}
		}
		return osst.str();
	}
};


//special case of the array2D constructor: if the pointer is a general arraydata<double> object
//useful to create 2D objects from general arraydata objects ...
template <> template <>
inline array2D<double>::array2D( arraydata<double> *data, unsigned int size_dim1, unsigned int size_dim2)
		: arraydata<double>( data )
{
	setDim( size_dim1, size_dim2 );
	if (data->size() != size_dim1*size_dim2) {
		std::cerr << "WARNING in array2D::array2D. Inconsistent array size. ";
		std::cerr << "size1D=" << data->size() << ", size2D=" << size_dim1*size_dim2 
			<< "=" << size_dim1 << "*" << size_dim2 << "" << std::endl;
	}
}

//... or from array1D objects
template <> template <>
inline array2D<double>::array2D( array1D<double> *data, unsigned int size_dim1, unsigned int size_dim2)
		: arraydata<double>( data )
{
	setDim( size_dim1, size_dim2 );
	if (data->size() != size_dim1*size_dim2) {
		std::cerr << "WARNING in array2D::array2D. Inconsistent array size. ";
		std::cerr << "size1D=" << data->size() << ", size2D=" << size_dim1*size_dim2 
			<< "=" << size_dim1 << "*" << size_dim2 << "" << std::endl;
	}
}



//=================================================================================
//
// array3D -- 3-dimensional array
//
//=================================================================================
template<typename T>
class array3D : public arraydata<T>{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;

	void setDim( unsigned int dim1, unsigned int dim2, unsigned int dim3 ){
		p_dim1 = dim1;
		p_dim2 = dim2;
		p_dim3 = dim3;
	}		
public:
	//-----------------------------------------------------------------------------------------constructors
	array3D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1, unsigned int size_dim3 = 1 )
			: arraydata<T>(size_dim1*size_dim2*size_dim3){
		setDim( size_dim1, size_dim2, size_dim3 );
	}

	//-----------------------------------------------------------------------------------------copy functions	
	void copy( const array3D& src ){
		setDim( src.dim1(), src.dim2(), src.dim3() );
		this->arraydata<T>::copy( src );
	}
	
	unsigned int dim1() const 						{ return p_dim1; }
	unsigned int dim2() const 						{ return p_dim2; }
	unsigned int dim3() const 						{ return p_dim3; }
	T get( unsigned int i, unsigned int j, unsigned int k ) const {
		assert(i < dim1());
		assert(j < dim2());
		assert(k < dim3());
		return arraydata<T>::get_atIndex( k*dim1()*dim2() + j*dim1() + i );
	}
	void set( unsigned int i, unsigned int j, unsigned int k, T value ){
		assert(i < dim1());
		assert(j < dim2());
		assert(k < dim3());
		arraydata<T>::set_atIndex( k*dim1()*dim2() + j*dim1() + i, value ); 
	}

    std::string getASCIIdata( bool annotate=1 ) const{
		std::ostringstream osst;
		if (this->size() == 0) {
			osst << "3D data has size zero." << std::endl;
		}else{
			if (annotate){ 
				osst << "3D data, dim1=" << dim1() << ", dim2=" << dim2() << ", dim3=" << dim3() 
					<< ", size=" << this->size() << std::endl; 
			}
			for (int k = 0; k<dim3(); k++){
				if (annotate){ osst << " [[" << std::endl; }
				for (int j = 0; j<dim2(); j++){
					if (annotate){ osst << "  ["; }
					for (int i = 0; i<dim1(); i++) {
						osst << " " << get(i, j, k);
					}//k
					if (annotate){ osst << "]"; }
					osst << std::endl;
				}//j
				if (annotate){ osst << "]]"; }
				osst << std::endl;
			}//i
		}
		return osst.str();
	}
};




//=================================================================================
//
// array4D -- 4-dimensional array
//
//=================================================================================
template<typename T>
class array4D : public arraydata<T>{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;
	int p_dim4;
	
	void setDim( unsigned int dim1, unsigned int dim2, unsigned int dim3, unsigned int dim4 ){
		p_dim1 = dim1;
		p_dim2 = dim2;
		p_dim3 = dim3;
		p_dim4 = dim4;
	}

public:
	//-----------------------------------------------------------------------------------------constructors
	array4D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1, 
				unsigned int size_dim3 = 1, unsigned int size_dim4 = 1 )
				: arraydata<T>(size_dim1*size_dim2*size_dim3*size_dim4){
		setDim( size_dim1, size_dim2, size_dim3, size_dim4 );
	}
	array4D( arraydata<T> *array, 
				unsigned int size_dim1, unsigned int size_dim2, 
				unsigned int size_dim3, unsigned int size_dim4 )
				: arraydata<T>(array){
		setDim( size_dim1, size_dim2, size_dim3, size_dim4 );
	}
    
	//-----------------------------------------------------------------------------------------copy functions
    void copy( const array4D<T>& src ){
		setDim( src.dim1(), src.dim2(), src.dim3(), src.dim4() );
		arraydata<T>::copy( src.data() );
	}
	
	unsigned int dim1() const			{ return p_dim1; }
	unsigned int dim2() const			{ return p_dim2; }
	unsigned int dim3() const			{ return p_dim3; }
	unsigned int dim4() const			{ return p_dim4; }
	
	double get( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const{
		assert(i < dim1());
		assert(j < dim2());
		assert(k < dim3());
		assert(l < dim4());
		return arraydata<T>::get_atIndex( l*dim1()*dim2()*dim3() + k*dim1()*dim2() + j*dim1() + i );
	}
	
	void set( unsigned int i, unsigned int j, unsigned int k, unsigned int l, double value ){
		assert(i < dim1());
		assert(j < dim2());
		assert(k < dim3());
		assert(l < dim4());
		arraydata<T>::set_atIndex( l*dim1()*dim2()*dim3() + k*dim1()*dim2() + j*dim1() + i, value);
	}

	void getRepresentationIn2D( array2D<double> *&img ){
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
		delete img;
		img = new array2D<double>(dim1()*dim4(), dim2()*dim3());	
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


	    std::string getASCIIdata( bool annotate=1 ) const{
		std::ostringstream osst;
		if (this->size() == 0) {
			osst << "4D data has size zero." << std::endl;
		}else{
			if (annotate){ 
				osst << "4D data, dim1=" << dim1() << ", dim2=" << dim2() 
					<< ", dim3=" << dim3() << ", dim4=" << dim4() << ", size=" << this->size() << std::endl; 
			}
			for (int l = 0; l<dim4(); l++){
				for (int k = 0; k<dim3(); k++){
					if (annotate){ osst << " [[" << std::endl; }
					for (int j = 0; j<dim2(); j++){
						if (annotate){ osst << "  ["; }
						for (int i = 0; i<dim1(); i++) {
							osst << " " << get(i, j, k, l);
						}//i
						if (annotate){ osst << "]"; }
						osst << std::endl;
					}//j
					if (annotate){ osst << "]]"; }
					osst << std::endl;
				}//k
			}//l
		}
		return osst.str();
	}
};


#endif
