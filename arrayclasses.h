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

#include <string>

//tell the compiler that these classes are going to be defined below (needed for copy constructors in arraydata)
class array1D;
class array2D;
class array3D;


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
    arraydata( const int16_t *CArray, const unsigned int size_val );     //initialize with C-array of ints of a given length
	arraydata( const float *CArray, const unsigned int size_val );			//initialize with C-array of floats of a given length
    arraydata( const arraydata &src );                          //copy constructor
    arraydata & operator=(const arraydata & src);               //assignment operator
	~arraydata();

    //helper functions
    void init();
    void copy( const double* src_data, const unsigned int size_val );
    void copy( const float* src_data, const unsigned int size_val );
    void copy( const int* src_data, const unsigned int size_val );
    void copy( const arraydata& src );
    void copy( const array1D& src );
    void destroy();
    double *data() const;                            			//return pointer to internal raw data array
    
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
	double calcMax() const;
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
};




//=================================================================================
//
// array1D
//
//=================================================================================
class array1D : public arraydata{
	
private:
	int p_dim1;
	
public:
	array1D( unsigned int size_dim1 = 1);                          
    array1D( int16_t *CArray, unsigned int size_val );          //init with C-style array of ints
	array1D( float *CArray, unsigned int size_val );			//init with C-style array of floats
    array1D( array2D* dataTwoD );                               //init with an array2D object
	~array1D();
    
    void copy( const array1D& src );
	
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	
	double get( unsigned int i ) const;
	void set( unsigned int i, double value );
    
    std::string getASCIIdata() const;
    int writeToASCII( std::string filename ) const;
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
	
public:
	array2D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1 );              //default constructor
    array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2);   // use 1D data to initialize
	~array2D();

    void copy( const array2D& src );
        
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2() const;
	void setDim2( unsigned int size_dim2 );
	
	double get( unsigned int i, unsigned int j ) const;                 //returns single pixel value
	void set( unsigned int i, unsigned int j, double value );

    std::string getASCIIdata() const;
	int writeToASCII( std::string filename, int format = 0 ) const;		// format default(0): 2D; (1): column output; (2): row output
	
	//-------------functions special to 2D case-----------------
    int getCol( int colnum, array1D *&col ) const;						//returns one dimensional column, extracted at the specified column number
	void setCol( int colnum, const array1D *col, int start=0 );			//sets a one-dimensional column, beginning at a 'start' value				
    
    int getRow( int rownum, array1D *&row ) const; 	                    //returns one-dimensional 'row' or 'col'
	void setRow( int rownum, const array1D *row, int start=0 );

	void transpose();													//transpose (dim1,dim2) --> (dim2,dim1)
	void flipud();														//flip up-down
	void fliplr();														//flip left-right
	
	void gradientAlongDim1( double lowlim, double highlim );			//linear gradient along one dimension, same along other dimension
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
	
public:
//    array3D();
	array3D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1, unsigned int size_dim3 = 1 );  //default constructor
	~array3D();
    
    void copy( const array3D& src );
	
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2() const;
	void setDim2( unsigned int size_dim2 );
	unsigned int dim3() const;
	void setDim3( unsigned int size_dim3 );
	
	double get( unsigned int i, unsigned int j, unsigned int k ) const;
	void set( unsigned int i, unsigned int j, unsigned int k, double value );

    std::string getASCIIdata() const;
    int writeToASCII( std::string filename ) const;
};


#endif
