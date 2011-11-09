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
#include <string>
#include <vector>

//tell the compiler that these classes are going to be defined below (needed for copy constructors in arraydata)
class array1D;
class array2D;
class array3D;
class array4D;


//=================================================================================
//
// arraydata (the base class for array1D, array2D, array3D)
//
//=================================================================================
class arraydata {
protected:
	double *p_data;												//pointer to the data array
	unsigned int p_size;	
	
public:
	arraydata( const unsigned int sizeval = 1 );
	arraydata( const arraydata *array );
	arraydata( const array1D *array );
	arraydata( const array2D *array );
	arraydata( const array3D *array );
	arraydata( const array4D *array );
	arraydata( const std::vector<double> vec );
	
	template <class T> arraydata( const T *CArray, const unsigned int size_val ){
		init();
		this->copy( CArray, size_val );
	}
	
	arraydata( const arraydata &src );                          //copy constructor
    arraydata & operator=(const arraydata & src);               //assignment operator
	~arraydata();

    //helper functions
    void init();												//initialize class
	void destroy();												//clean up class memory
	double *data() const;                            			//return pointer to internal raw data array
    
	//functions to create copies from various sources
	void copy( const arraydata *src );
	void copy( const array1D *src );
	void copy( const array2D *src );
	void copy( const array3D *src );
	void copy( const array4D *src );
	void copy( const std::vector<double> vec );
	
	void copy( const arraydata& src );
	
	// templated function for c-style arrays of basic types (float, double, int, ...)
	template <class T> void copy( const T *src, const unsigned int arraysize ){
		p_size = arraysize;
		if (p_size > 0 && src){
			this->destroy();
			p_data = new double[ p_size ];
			for (unsigned int i = 0; i < p_size; i++) {
				p_data[i] = (double)src[i];
			}
		}else{
			p_data = NULL;
			std::cerr << "Error in arraydata::copy! src=" << src << ", size=" << arraysize << std::endl;
		}
	}
    
    void zeros();												//set all elements to 0
    void zeros( unsigned int start, unsigned int stop );		//set all elements between start and stop to 0 (incl. start, excl. stop)
	void ones();												//set all elements to 1
    void ones( unsigned int start, unsigned int stop );			//set all elements between start and stop to 1 (incl. start, excl. stop)
	void range( double neg, double pos );						//set elements to a range of values, given by the boundaries
    
	// 'atIndex' functions:
	// the following functions change or return properties
	// assuming a one-dimensional array (consistent with the internal data structure), 
    // no matter what dimension the actual subclass may have
	unsigned int size() const;										//total size of the array	
	double get_atIndex( unsigned int index) const;					// get element value
	void set_atIndex( unsigned int index, double val);			// set element value
	
	double calcMin() const;
	double calcMin(int &pos) const;
	double calcMax() const;
	double calcMax(int &pos) const;
	double calcSum() const;
	double calcAvg() const;
	
    std::string getASCIIdataAsRow() const;                         //can/should be overridden by subclasses
    std::string getASCIIdataAsColumn() const;
    
    //perform some basic math on array
    int addValue( double val );
	int subtractValue( double val );
    int multiplyByValue( double value );
	int divideByValue( double value );
 
	int addArrayElementwise( const arraydata *secondArray );
	int subtractArrayElementwise( const arraydata *secondArray );
    int multiplyByArrayElementwise( const arraydata *secondArray );
    int divideByArrayElementwise( const arraydata *secondArray );
	
	// application of a mask (usually ones and zeros)
	// set values below 'checkval' to 'rejectval'
	int applyMask( arraydata* mask, double checkval = 0.5, double rejectval = 0. );	
	
	int getHistogram( array1D *&hist, array1D *&bins, unsigned int nBins = 20 );
	std::string getHistogramASCII( unsigned int nBins = 20 );

	int getHistogramInBoundaries( array1D *&hist, array1D *&bins, unsigned int nBins, double min, double max );
	std::string getHistogramInBoundariesASCII( unsigned int nBins, double min, double max );
};




//=================================================================================
//
// array1D
//
//=================================================================================
class array1D : public arraydata{
	
private:
	int p_dim1;

	//private setters for dimensions
	//not (yet) safe to resize this class by changing dimensions after instatiation
	void setDim1( unsigned int size_dim1 );	
			
public:
	array1D( unsigned int size_dim1 = 1);                          
    array1D( arraydata* data );                               	//init with any arraydata object
	array1D( array1D* dataOneD );								//init with an array1D object
    array1D( array2D* dataTwoD );                               //init with an array2D object
	template <class T> array1D( T *CArray, unsigned int size_dim1 )
        : arraydata( CArray, size_dim1 ){
		setDim1( size_dim1 );
	}
		
	~array1D();
    
    void copy( const array1D& src );
	
	unsigned int dim1() const;
	
	double get( unsigned int i ) const;
	void set( unsigned int i, double value );
    
    std::string getASCIIdata( bool annotate=1 ) const;

};




//=================================================================================
//
// array2D
//
// this class, while as generic as possible, adopts the convention (rows, columns)
// of packages like matlab or python for its arguments
//=================================================================================
class array2D : public arraydata{
private:
	int p_dim1;
	int p_dim2;

	//private setters for dimensions
	//not (yet) safe to resize this class by changing dimensions after instatiation
	void setDim1( unsigned int size_dim1 );	
	void setDim2( unsigned int size_dim2 );	
			
public:
	array2D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1 );              //default constructor
    array2D( arraydata* data, unsigned int size_dim1, unsigned int size_dim2);
	array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2);	// use 1D data to initialize
	array2D( array2D* dataTwoD );													// initialize with a copy of the argument
	template <class T> array2D( T *dataCArray, unsigned int size_dim1, unsigned int size_dim2 )
			: arraydata( dataCArray, size_dim1*size_dim2 ){
		setDim1( size_dim1 );
		setDim2( size_dim2 );	
	}

	~array2D();

    void copy( const array2D& src );
        
	unsigned int dim1() const;
	unsigned int dim2() const;
	
	double get( unsigned int i, unsigned int j ) const;                 //returns single pixel value
	void set( unsigned int i, unsigned int j, double value );

    std::string getASCIIdata( bool annotate=1 ) const;
	
	//-------------functions special to 2D case-----------------
    int getCol( int colnum, array1D *&col ) const;						//returns one dimensional column, extracted at the specified column number
	void setCol( int colnum, const array1D *col, int start=0 );			//sets a one-dimensional column, beginning at a 'start' value				
    
    int getRow( int rownum, array1D *&row ) const; 	                    //returns one-dimensional 'row' or 'col'
	void setRow( int rownum, const array1D *row, int start=0 );
	
	int calcAvgRow( array1D *&row ) const;
	int calcAvgCol( array1D *&col ) const;
	
	void transpose();													//transpose (dim1,dim2) --> (dim2,dim1)
	void flipud();														//flip up-down
	void fliplr();														//flip left-right
	
	
	//linear gradient along one dimension, same along other dimension
	void gradientAlongDim1( double lowlim, double highlim );			
	void gradientAlongDim2( double lowlim, double highlim );
	
	void generateTestPattern( int type );
};




//=================================================================================
//
// array3D
//
//=================================================================================
class array3D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;

	//private setters for dimensions
	//not (yet) safe to resize this class by changing dimensions after instatiation
	void setDim1( unsigned int size_dim1 );	
	void setDim2( unsigned int size_dim2 );	
	void setDim3( unsigned int size_dim3 );	
		
public:
	array3D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1, unsigned int size_dim3 = 1 );  //default constructor
	~array3D();
    
    void copy( const array3D& src );
	
	unsigned int dim1() const;
	unsigned int dim2() const;
	unsigned int dim3() const;
	
	double get( unsigned int i, unsigned int j, unsigned int k ) const;
	void set( unsigned int i, unsigned int j, unsigned int k, double value );

    std::string getASCIIdata( bool annotate=1 ) const;
};




//=================================================================================
//
// array4D
//
//=================================================================================
class array4D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;
	int p_dim4;
	

	//private setters for dimensions
	//not (yet) safe to resize this class by changing dimensions after instatiation
	void setDim1( unsigned int size_dim1 );	
	void setDim2( unsigned int size_dim2 );	
	void setDim3( unsigned int size_dim3 );	
	void setDim4( unsigned int size_dim4 );

public:
	array4D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1, 
				unsigned int size_dim3 = 1, unsigned int size_dim4 = 1 );  //default constructor
	array4D( array1D *dataOneD, 
				unsigned int size_dim1, unsigned int size_dim2, 
				unsigned int size_dim3, unsigned int size_dim4 );
	~array4D();
    
    void copy( const array4D& src );
	
	unsigned int dim1() const;
	unsigned int dim2() const;
	unsigned int dim3() const;
	unsigned int dim4() const;
	
	double get( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const;
	void set( unsigned int i, unsigned int j, unsigned int k, unsigned int l, double value );

	void getRepresentationIn2D( array2D *&img );

    std::string getASCIIdata( bool annotate=1 ) const;
};


#endif
