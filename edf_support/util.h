//========================================================================
// File: util.h
//
// Author: Christian G. Schroer
// Description: This file contains class descriptions for numerous 
//		statistical functions. 
//========================================================================

#ifndef util_h_flag
#define util_h_flag

#include <vector>

#include <string>
//#include <sstream>
//#include <cmath>
//#include <cstdlib>
//#include <cstdio>
//#include <time.h>

#include <iostream>

#include "timer.h"

namespace matlibutil
{
	//--------------------------------------------------------------------------
	// SIGMA is the standard deviation of the y-coordinate in the linear fit.
	// Right now, all are assumed to be the same
	// RAND_MAX = 2147483647, defined in cstdlib
	//--------------------------------------------------------------------------
	const double SIGMA = 1; 

	#ifdef WIN32
		typedef unsigned int useconds_t; // for usleep()
		
		//typedef unsigned int DWORD;	
			//by JF 09/2008: removed this line, because DWORD is already defined in WinDef.h
			//
	#endif

	//==========================================================================
	//                                                             Stream output
	//==========================================================================

	//--------------------------------------------------------------------------
	// This defines a new class C_Binary for stream output formatting.
	// The intention is to have a possibility to print numeric values
	// in binary format (oct, dec and hex formats are already defined
	// in standard library).
	// Expample:
	// unsigned int value = 0x12345678;
	// cout << std::hex << value; --> 12345678                        (built-in)
	// cout << std::dec << value; --> 305419896                       (built-in)
	// cout << std::oct << value; --> 2215053170                      (built-in)
	// cout << matlibutil::bin(value);   --> 00010010001101000101011001111000 (new)
	//--------------------------------------------------------------------------
	template <class T>
	class C_Binary
	{
		public:
			C_Binary()
			{
				p_value=0;
			}

			C_Binary( T value )
			{
				p_value = value;
			}

			inline T highbit(T &value)
			{
				return value = (((T)(-1)) >> 1) + 1;
			}

			std::ostream &toBin( std::ostream &os, T& value )
			{
				for ( T bit = highbit(bit); bit; bit >>= 1 )
				{
					os << ( ( value & bit ) ? '1' : '0' );
				}
				return os;
			}

			std::ostream &toBin( std::ostream &os )
			{
				for ( T bit = highbit(bit); bit; bit >>= 1 )
				{
					os << ( ( p_value & bit ) ? '1' : '0' );
				}
				return os;
			}

		private:
			T p_value;
	};

	//--------------------------------------------------------------------------
	// Output stream operator "<<" for our C_Binary class.
	// Example:
	// std::C_Binary<unsigned int> binary(0x12345678);
	// cout << binary;    ---> 00010010001101000101011001111000
	//--------------------------------------------------------------------------
	template <class T>
	std::ostream &operator<<( std::ostream &os, C_Binary<T> binary )
	{
		return binary.toBin( os );
	}

	//--------------------------------------------------------------------------
	// Function that returns an instance of Class C_Binary
	// with appropriate type.
	// Example:
	// bin(0x12345678);   ---> instance of C_Binary with p_value==0x12345678.
	//--------------------------------------------------------------------------
	template <class T>
	C_Binary<T> bin( T value )
	{
		return C_Binary<T>(value);
	}


	//==========================================================================
	//                                                               Simple math
	//==========================================================================

	//--------------------------------------------------------------------------
	// Calculates the representation of <value> with respect to place-value
	// system of given base. The representation is stored within <digits>
	// The first integer (digits[0]) will be the place with the lowest value,
	// the last integer (digits[digits.size()-1) will be the place with the
	// highest value, and it will be zero if and only if <value> is zero.
	// <digits> will be empty if <base> is smaller than two.
	// The sign of value will be ignored.
	//--------------------------------------------------------------------------
	void placeValueSystem( double value, int base, std::vector<int> &digits );

	//--------------------------------------------------------------------------
	// Calculates the representation of <value> with respect to place-value
	// system of given base.
	// The integral part of the representation (the digits to the left of
	// the fractional point) is stored within <integral_digits>.
	// The fractional part of the representation (the digits to the right of
	// the fractional point) is stored within <fractional_digits>.
	// The first element of <integral_digits> will be the place with the
	// lowest value, the last element of <integral_digits> will be the place
	// with the highest value, and it will be zero if and only if <value> is
	// in the open interval (-1, 1) = [-1, 1]\{-1, 1}.
	// <integral_digits> will be empty (zero sized) if and only if <base> is
	// smaller than 2.
	// The n-th element of <fractional_digits> is the n-th place to the left
	// of the fractional point. The size of <fractional_digits> is given by
	// <precision>. The direction of rounding is given by the nearest distance
	// of one of the two possible directions. If the distance is zero, no
	// rounding will be performed, if the distance is equal, rouning will be
	// done to that number the absolute of which is greater.
	// <fractional_digits> will be empty (zero sized), if <precision> is
	// smaller than 1, or if <base> is smaller than 2.
	// The sign of value will be ignored.
	//--------------------------------------------------------------------------
	void placeValueSystem( double value, int base, int precision,
						   std::vector<int> &integral_digits,
						   std::vector<int> &fractional_digits );

	//--------------------------------------------------------------------------
	// Calculates the representation of <value> with respect to place-value
	// system of given base and returns the digit at given place.
	// Places of the integral part of the representation (the digits to the left
	// of the fractional point) are counted with positive numbers (starting with
	// zero for the less significant place) whereas places of the fract-
	// ional part of the representation are counted with negative number.
	// If <precision> is greater than zero, the rounding of fractional digits
	// will be performed with that precision. Otherwise, rounding will be done
	// with a precision dependent on the value of <place>. If place is negative
	// (fractional place), precision will be given by -place+1. Otherwise,
	// precision will be set to 1.
	// Examples:
	// getDigit( 56789.01234, 10,  0 ) --> 9
	// getDigit( 56789.01234, 10,  4 ) --> 5
	// getDigit( 56789.01234, 10, -1 ) --> 0
	// getDigit( 56789.01234, 10, -5 ) --> 4
	// getDigit(    21      ,  2,  0 ) --> 1
	// getDigit(    21      ,  2,  1 ) --> 0
	// getDigit(    21      ,  2, -5 ) --> 0
	//--------------------------------------------------------------------------
	int getDigit( double value, int base, int place, int precision=0 );

	//--------------------------------------------------------------------------
	// Searches for a solution of f(x)=y using Newton iteration.
	// Starts at xstart and iterates as long as |f(x)-y|>|delta|.
	// g is the derivative of f.
	// The vector 'parameters' defines optional parameters.
	// Does only work if xstart is within the convergence interval.
	// (Should be improved with respect to this matter.)
	//--------------------------------------------------------------------------
	typedef double (*newtonfunc_t)(double x, const std::vector<double> &parameters);
	double newtonIteration( newtonfunc_t f, newtonfunc_t g,
							const std::vector<double> &parameters,
							double y, double xstart, double delta );
	
	//--------------------------------------------------------------------------
	// Calculates the modulo of i with respect to n.
	// In contrast to the built-in operator '%', pmod will always have positive
	// results. While '-2%5' will result in '-2', 'pmod(-2,5)' will return '3'.
	//--------------------------------------------------------------------------
	int pmod( int i, int n );


	//==========================================================================
	//                                                                Byte order
	//==========================================================================

	int BigEndianByteOrder(void);

	//--------------------------------------------------------- reverseByteOrder
	// Reverses the byte order of the integer numbers in a data field.
	// char *data  : data field to be modified
	// sizeDataType: size of the datatype (in bytes)
	// length      : length of the data field
	//--------------------------------------------------------------------------
	void reverseByteOrder( char *data, int sizeDataType, int length=1 );

	void revByteOrder( float			*fval );
	void revByteOrder( double			*dval );
	void revByteOrder( short			*sval );
	void revByteOrder( int				*ival );
	void revByteOrder( long 			*lval );
	void revByteOrder( unsigned short	*sval );
	void revByteOrder( unsigned int 	*ival );
	void revByteOrder( unsigned long	*lval );
	void revByteOrder( float			*data, unsigned long length );
	void revByteOrder( double			*data, unsigned long length );
	void revByteOrder( short			*data, unsigned long length );
	void revByteOrder( int				*data, unsigned long length );
	void revByteOrder( long 			*data, unsigned long length );
	void revByteOrder( unsigned short	*data, unsigned long length );
	void revByteOrder( unsigned int 	*data, unsigned long length );
	void revByteOrder( unsigned long	*data, unsigned long length );
	unsigned short convlitend(unsigned short &var);

	//==========================================================================
	//                                                         vector operations
	//==========================================================================

	double scalarProd(std::vector<double> &v1, std::vector<double> &v2);

	//==========================================================================
	//                                                      random distributions
	//==========================================================================

	double equipart(double low, double high);
	int equipart(int low, int high);
	unsigned int equipart(unsigned int low, unsigned int high);
	int exprand(double E);
	//long int poisson(double average);
	int randsphere(double& x, double& y, double& z);
	double poisson(double average);
	void boxmueller(double& x1, double& x2);
	void gauss(double& x1, double& x2, double average, double stddev);
	double gauss(double average, double stddev);
	double gausstrunc(double average, double stddev, double highlimit, double lowlimit);
	double gauss(double average);
	int selectBin(std::vector<double> width);
	int bisect(std::vector<double> array);
	void heapsort(std::vector<double> &array);
	void indextable(std::vector<double> array, std::vector<unsigned int> &indexVec);
	void init_rand();
	void isotrope(double& theta, double& phi);


	void complexSqrt(double &re, double &im);
	double erfcc(double x);
	double erfkn(double x);
	double fermiFct(double x, double prec);



	//==========================================================================
	//                                                         string formatting
	//==========================================================================

	std::string appendNumber(std::string base, int number, int number_digits = 4);
	std::string appendNumber(std::string base, unsigned int number, int number_digits = 4);
	std::string appendNumber(std::string base, long int number, int number_digits = 4);
	std::string appendNumber(std::string base, long unsigned int number, int number_digits = 4);
	std::string appendNumber(std::string base, double number, int length = 6, int num_digits = 5);
	std::string appendNumber(std::string left, std::string number, std::string right, std::string extension, int digits = 4);
	int appendNumber(std::string& base, std::string left, int number, std::string right, int number_digits);

	std::string toString(std::string number);
	std::string toString(int number);
	std::string toString(int number, int digits);
	std::string toString(long int number);
	std::string toString(unsigned int number);
	std::string toString(unsigned long int number);
	std::string toString(double number, int prec = 8);
	int toInt(std::string number);
	double toDouble(std::string number);

	std::string toLower( std::string text );
	std::string toUpper( std::string text );

	//--------------------------------------------------------------------------
	// Compares the strings A and B case insensitively.
	// Returns true, if A and B are equal apart from lower/upper case.
	//--------------------------------------------------------------------------
	bool equal( std::string A, std::string B );

	//--------------------------------------------------------------------------
	// Returns a string, so that 'text+string' has the length given by width.
	// Returns the empty string, if text is longer than width.
	// As an alternative, the filling character can be changed to any ascii
	// character. This can be done with <fill>.
	//--------------------------------------------------------------------------
	std::string spacepadding( std::string text, int width, char fill=' ' );

	void set_keypress();
	void reset_keypress();

	int get_cursorPosX();
	int get_cursorPosY();
	void get_cursorPos( int &row, int &col );

	std::string clearLine();
	std::string clearLineToEnd();
	std::string clearLineToBegin();
	std::string clearScreen();
	std::string clearScreenToEnd();
	std::string clearScreenToBegin();
	std::string cursorUp( int num_lines=1 );
	std::string cursorDown( int num_lines=1 );
	std::string cursorLeft( int num_cols=1 );
	std::string cursorRight( int num_cols=1 );
	std::string cursorToPosition( int row, int col );
	std::string cursorToCol( int col );
	std::string cursorToRow( int row );
	std::string cursorLineBegin();
	std::string cursorLineEnd();
	std::string hideCursor();
	std::string showCursor();
	std::string textreset();
	std::string boldon();
	std::string boldoff();
	std::string fainton();
	std::string faintoff();
	std::string italicon();
	std::string italicoff();
	std::string underline1on();
	std::string underline2on();
	std::string underlineoff();
	std::string blinkslowon();
	std::string blinkfaston();
	std::string blinkoff();
	std::string inverton();
	std::string invertoff();
	std::string hideon();
	std::string hideoff();
	std::string strikethruon();
	std::string strikethruoff();

	std::string ANSIescape( std::string code );


	//==========================================================================
	//                                                    file system management
	//==========================================================================
	
	//--------------------------------------------------------------------------
	// Returns platform dependant path tree separator ( '/' or '\' ).
	//--------------------------------------------------------------------------
	std::string get_separator();
	
	//--------------------------------------------------------------------------
	// Needed to convert path names from linux to windows convention.
	//--------------------------------------------------------------------------
	std::string slashesToBackslashes( std::string path );

	//--------------------------------------------------------------------------
	// Needed to convert path names from windows to linux convention.
	//--------------------------------------------------------------------------
	std::string backslashesToSlashes( std::string path );

	//--------------------------------------------------------------------------
	// Unix/Apple: convert backslashes to slashes.
	// MS Windows: convert slashes to backslashes.
	//--------------------------------------------------------------------------
	std::string properPath( std::string path );

	//--------------------------------------------------------------------------
	// Receive a path and return only the directory part of it.
	// Works with slashes or backslashes.
	// Example: "this/is/it" -> "this/is/"
	// Example: "this\is\it" -> "this\is\"
	//--------------------------------------------------------------------------
	std::string returnDirectoryPartOf( std::string path, bool verbose=false);

	//--------------------------------------------------------------------------
	// Receive a path and return only the filename part of it.
	// Works with slashes or backslashes,
	// Example: this/is/it -> it
	// Example: this\is\it -> it
	//--------------------------------------------------------------------------
	std::string returnFilenamePartOf( std::string path, bool verbose=false );

	//--------------------------------------------------------------------------
	// Extracts the parts from a path tree and returns them within a vector of
	// strings. An absolute path results in a vector with the empty string as
	// the first element. Consecutive separators are treated as just one separator.
	// Example: "aaa/bbb/ccc/ddd" ---> ( "aaa", "bbb", "ccc", "ddd" )
	// Example: "/aaa/bbb/ccc"    ---> ( "", "aaa", "bbb", "ccc" )
	// Example: "aaa///bbb/ccc/"  ---> ( "aaa", "bbb", "ccc" )
	//--------------------------------------------------------------------------
	int resolvePathTree( std::string pathtree, std::vector<std::string> &parts );
	int resolvePathTree( std::string pathtree, std::vector<std::string> &parts, char separator );

	//--------------------------------------------------------------------------
	// Creates the full path tree, even if more than the last part of the path
	// does not exist. Fails, if any of the parts of the path tree already
	// exists as a non-directory file type.
	// Example: 'aaa/bbb' + createPathTree( "aaa/bbb/ccc/ddd" ) --> 'aaa/bbb/ccc/ddd'
	// Example: 'aaa/bbb' + createPathTree( "aaa" ) --> 'aaa/bbb'
	//--------------------------------------------------------------------------
	int createPathTree( std::string pathtree, bool verbose=false );

	//--------------------------------------------------------------------------
	// Returns true if 'fullname' does exist as a regular file.
	// Returns false if 'fullname' does not exist or if 'fullname'
	// is not a regular file (e.g. a directory or a fifo).
	//--------------------------------------------------------------------------
	bool fileDoesExist( std::string fullname );

	//--------------------------------------------------------------------------
	// Returns true if 'fullname' does exist as a fifo file.
	// Returns false if 'fullname' does not exist or if 'fullname'
	// is not a fifo file (e.g. a directory or a regular file).
	//--------------------------------------------------------------------------
	bool fifoDoesExist( std::string fullname );

	//--------------------------------------------------------------------------
	// Returns true if 'path' does exist as a directory.
	// Returns false if 'fullname' does not exist or if 'fullname'.
	// is not a directory.
	//--------------------------------------------------------------------------
	bool directoryDoesExist( std::string path, bool verbose=false );

	//--------------------------------------------------------------------------
	// Checks, if the file exists in the current directory.
	// Returns true if file exists and can be opened.
	// In contrast to 'fileDoesExist()', this functions fails,
	// if the given file exists but cannot be opened because
	// of missing permissions or something similar.
	//--------------------------------------------------------------------------
	bool checkFileExistance(std::string FileName);

	//--------------------------------------------------------------------------
	// Returns the number of lines in a file.
	//--------------------------------------------------------------------------
	long int countLines(std::string fn);
	
	//--------------------------------------------------------------------------
	// Returns the number of lines in a file.
	//--------------------------------------------------------------------------
	long int countLines(std::string fn, std::string skipchars);



	//==========================================================================
	//                                                               portability
	//==========================================================================

	//--------------------------------------------------------------------------
	// Platform independent implementation of the posix usleep() function.
	// On unix:  calls the built-in usleep() declared in unistd.h.
	// On WIN32: calls the built-in Sleep() with a ms<-->us conversion.
	//--------------------------------------------------------------------------
	int usleep( useconds_t time_usec );



	//==========================================================================
	//                                                                statistics
	//==========================================================================

	class binio
	{
		public:
			binio();
			~binio();
			int write(FILE *fp, char *buffer, int size);
			int read(FILE *fp, char *buffer, int size);
	};

	class Min
	{
		public:
			Min();
			~Min();
			void reset();
			double get();
			int compare(double x);
		private:
			double minimum;
	};

	class Max
	{
		public:
			Max();
			~Max();
			void reset();
			double get();
			int compare(double x);
		private:
			double maximum;
	};
	
	class Range
	{
		public:
			Range();
			~Range();
			void reset();
			double get();
			double getMax();
			double getMin();
			int compare(double x);
		private:
			double maximum, minimum;
	};

	class IMin
	{
		public:
			IMin();
			~IMin();
			void reset();
			long int get();
			int compare(long int x);
		private:
			long int minimum;
	};

	class IMax
	{
		public:
			IMax();
			~IMax();
			void reset();
			long int get();
			int compare(long int x);
		private:
			long int maximum;
	};

	class UIMax
	{
		public:
			UIMax();
			~UIMax();
			void reset();
			unsigned long int get();
			int compare(unsigned long int x);
		private:
			unsigned long int maximum;
	};
	
	class sum
	{
		public:
			sum();
			~sum();
			void reset();
			double get();
			long int getCount();
			void add(double x);
		private:
			double summe;
			long int count;
	};

	class mean
	{
		public:
			mean();
			~mean();
			void reset();
			double get();
			double getSum();
			long int getCount();
			void add(double x, unsigned int n = 1);
		private:
			double sum;
			long int count;
	};
	
	
	class histogram
	{
		public:
			histogram(double lowval, double highval, int intervals);
			histogram(histogram *histoval);
			void add(histogram *histoval);
			void add(double val);
			void add(unsigned int val);
			void add_nocheck(unsigned int val);
			void clear();
			unsigned long int get(unsigned int bin);
			unsigned long int get(double val);
			void get(std::vector<unsigned long int > &histval);
			void get(std::vector<double > &histval);
			double getHigh();
			double getLow();
			unsigned long int getSize();
			double getMean(int center = 1); // center == 1: center of bin, center = 0: lower bound
			double getStDev();
			double getErrorOfMean();// Mittlerer Fehler des Mittelwerts
			double getPercentile(double percentval);
			double getMedian();
			double getErrorOfMedian();// Standardabweichung des Medians
			double getErrorOfMedian(unsigned int trials);
			void setHigh(double highval);
			void setLow(double lowval);
			unsigned long int getCount();
			unsigned long int sum();
		private:
			std::vector<unsigned long int > hist;
			double low;
			double high;
			unsigned long int count;
	};

	class qMoment
	{
		public:
			qMoment();
			~qMoment();
			void setMoment(long int countval, double momentval);
			void addMoment(long int countval, double momentval);
			int setQ(double qval);
			void reset();
			double get();
			double getQPow();
			long int getCount();
			void add(double x);
		private:
			double sum;
			double q;
			long int count;
	};

	class counter
	{
		public:
			counter();
			~counter();
			void reset();
			long int getCount();
			void add();
		private:
			long int count;
	};

	class Flag
	{
		public:
			Flag(double comparatorval);
			~Flag();
			void reset();
			void setComparator(double comparatorval);
			double getComparator();
			void set(int flagval);
			int get();
			double compare(double x);
		private:
			int flag;
			double comparator;
	};

	class Stat
	{
		public:
			Stat();
			~Stat();
			void set(Stat *statval); 
			void add(Stat *statval);
			void reset();
			double getStd();
			double getStd2();
			double getMean();
			long int getCount();
			void add(double x);
			double getQuad();
			double getSum();
		private:
			double sum;
			double quad;
			long int count;
	};

	class statf
	{
		public:
			statf();
			~statf();
			void set(statf *statval);
			void add(statf *statval);
			void reset();
			float getStd();
			float getStd2();
			float getMean();
			long int getCount();
			void add(double x);
			float getQuad();
			float getSum();
		private:
			float sum;
			float quad;
			long int count;
	};

	class twodstat
	{
		public:
			twodstat();
			~twodstat();
			void set(twodstat *statval);
			void add(twodstat *statval);
			void reset();
			double getStd();
			double getStd2();
			double getMeanX();
			double getMeanY();
			long int getCount();
			void add(double x, double y);
			double getSum1();
			double getSum2();
			double getQuad();
		private:
			double sum1;
			double sum2;
			double quad;
			long int count;
	};

	class twodstatf
	{
		public:
			twodstatf();
			~twodstatf();
			void set(twodstatf *statval);
			void add(twodstatf *statval);
			void reset();
			float getStd();
			float getStd2();
			float getMeanX();
			float getMeanY();
			long int getCount();
			void add(float x, float y);
			float getSum1();
			float getSum2();
			float getQuad();
			float sum1;
			float sum2;
			float quad;
			long int count;
	};

	class stDev
	{
		public:
			stDev();
			~stDev();
			void set(stDev *devval);
			void add(stDev *devval);
			void reset();
			double get();
			double getErrorOfMean();// Mittlerer Fehler des Mittelwerts
			long int getCount();
			void add(double x, unsigned int n = 1);
			double getQuad();
			double getSum();
		private:
			double sum;
			double quad;
			long int count;
	};

	class C_Correlator
	{
		public:
			C_Correlator();
			C_Correlator(const C_Correlator &);
			~C_Correlator();
			
			C_Correlator& operator=(const C_Correlator &);
			void set(C_Correlator *devval);
			void add(C_Correlator *devval);
			void reset();
			double get();
			long int getCount();
			void add(double x1, double x2);
			double getQuad();
			double getSum();

		private:
			double sum;
			double quad;
			long int count;
	};

	class stdev2
	{
		public:
			stdev2();
			~stdev2();
			void set(stdev2 *devval);
			void add(stdev2 *devval);
			void reset();
			double get();
			long int getCount();
			void add(double x);
			double getQuad();
			double getSum();
		private:
			double sum;
			double quad;
			long int count;
	};

	class statS
	{
		public:
			statS();
			~statS();
			int reset();
			long int getCount();
			int add();
			int add(double sig);
			double get();
		private:
			double S;
			long int count;
	};

	class statSx
	{
		public:
			statSx();
			~statSx();
			int reset();
			long int getCount();
			int add(double x);
			int add(double x, double sig);
			double get();
		private:
			double S;
			long int count;
	};

	class statSy
	{
		public:
			statSy();
			~statSy();
			int reset();
			long int getCount();
			int add(double y);
			int add(double y, double sig);
			double get ();
		private:
			double S;
			long int count;
	};

	class statSxx
	{
		public:
			statSxx();
			~statSxx();
			int reset();
			long int getCount();
			int add(double x);
			int add(double x, double sig);
			double get();
		private:
			double S;
			long int count;
	};

	class statSxy
	{
		public:
			statSxy();
			~statSxy();
			int reset();
			long int getCount();
			int add(double x, double y);
			int add(double x, double y, double sig);
			double get();
		private:
			double S;
			long int count;
	};

	class linFit
	{
		public:
			linFit();
			int add(double xval, double yval);
			int add(double xval, double yval, double sig);
			double getA();
			double getB();
			double getSlope();
			double getOffset();
		private:
			statS S;
			statSx Sx;
			statSy Sy;
			statSxx Sxx;
			statSxy Sxy;
			double Delta;
			double a;
			double b;
	};

	class Queue
	{
		public:
			Queue(long int lengthval);
			Queue(long int lengthval, double initval);
			~Queue();
			long int getLength();
			long int getCount();
			double getVal();
			double getVal(long int i);
			void add(double val);
			double getMean();
			double getDiff();
			double getLast();
		private:
			long int length;
			double *value;
			long int counter; // number of items sent to queue
	};


	
	//==========================================================================
	//                                                             interpolation
	//==========================================================================
	
	double interpolate(double val0, double val1, double t);
	double interpolate2D(double val00, double val01, double val10, double val11, 
						 double ty, double tx);
	
	// utils to convert from re/im to abs/phase and vice versa
	void absphase2reim(double abs, double phase, double &re, double &im);
	void reim2absphase(double re, double im, double &abs, double &phase);

	// interpolate complex numbers as re/im and abs/phase, respectively
	void interpolate(double valre0, double valim0, double valre1, double valim1, 
					 double t, double &resre, double &resim);
	void interpolatePhase(double valre0, double valim0, double valre1, double valim1, 
						  double t, double &resre, double &resim);

	// interpolate complex numbers in 2D, both in re/im and in abs/phase, respectively
	void interpolate2D(double valre00, double valim00, double valre01, double valim01, 
					   double valre10, double valim10, double valre11, double valim11, 
					   double ty, double tx, double &resre, double &resim);
	void interpolate2DPhase(double valre00, double valim00, double valre01, double valim01, 
							double valre10, double valim10, double valre11, double valim11, 
							double ty, double tx, double &resre, double &resim);


	//==========================================================================
	//                                                             miscellaneous
	//==========================================================================

//***commented out the following class for the use in the cross correlator project, but matlib generally needs this
/*
	// -------------------------------------------------------------- C_Progress
	//
	const double DEFAULT_TIMESTEP = 500; // in msec
	const char DEFAULT_PREFIX[] = "Progress: ";

	class C_Progress
	{
		public:
			C_Progress();
			C_Progress( double max_value );
			C_Progress( std::istream &in );
			C_Progress( std::ostream &out );
			C_Progress( const C_Progress &that );
			~C_Progress();

			C_Progress &operator=( const C_Progress &that );

			void set_timestep( double timestep );
			void set_prefix( std::string prefix );

			double get_timestep();
			std::string get_prefix();

			void init();
			void init( double max_value );
			void init( std::istream &in );
			void init( std::ostream &out );

			void show( int units_per_asterisk, int asterisks_per_line );
			void show( double value );
			void show( std::istream &in );
			void show( std::ostream &out );
		
			void finished();
		private:
			double				p_progress_0;
			double				p_progress_1;
			int 				p_progress;
			double				p_timestep;
			ns_timer::C_Timer	p_timer;
			std::string 		p_prefix;
	};
*/	
	
	// -------------------------------------------------------- C_LoopDataSingle
	//
	typedef double (*loopfunc_t)(double);
	
	class C_LoopDataSingle
	{
		public:
			C_LoopDataSingle();
			C_LoopDataSingle( double start, double stop, int steps );
			C_LoopDataSingle( double start, double stop, int steps, loopfunc_t func );
			C_LoopDataSingle( const C_LoopDataSingle &loopdatasingle );
			~C_LoopDataSingle();

			C_LoopDataSingle &operator=( const C_LoopDataSingle &loopdatasingle );
		
			void init();
		
			void set_start( double start );
			void set_stop( double stop );
			void set_steps( int steps );
			void set_func( loopfunc_t func );
			void set( double start, double stop, int steps );
			void set( double start, double stop, int steps, loopfunc_t func );
			void reset();
			
			double get_start() const;
			double get_stop() const;
			int get_steps() const;
			loopfunc_t get_func() const;
			
			bool next();
			double get_value() const;
			double get_sum( int depth=0 ) const;
			double get_sum_range( int start=0, int stop=-1 );
			int get_count() const;
			
		private:
			double p_start;
			double p_stop;
			int p_steps;
			int p_count;
			loopfunc_t p_func;
			std::vector<double> p_values;
	};


	// -------------------------------------------------------------- C_LoopData
	//
	typedef std::vector<C_LoopDataSingle> looplist_t;

	class C_LoopData
	{
		public:
			C_LoopData();
			C_LoopData( const C_LoopData &loopdata );
			~C_LoopData();
			
			C_LoopData &operator=( const C_LoopData &loopdata );
		
			void init();
			
			void clear();
			void reverse();
			void add( const C_LoopDataSingle &loopdatasingle );
			void add( double start, double stop, int steps );
			void add( double start, double stop, int steps, loopfunc_t func );
			void add_loops( looplist_t &looplist );
			
			C_LoopDataSingle get_loop_copy( int index ) const;
			C_LoopDataSingle &get_loop( int index );
			
			double get_start( int index ) const;
			double get_stop( int index ) const;
			int get_steps( int index ) const;
			loopfunc_t get_func( int index ) const;
			
			double get_value( int index ) const;
			double get_sum( int index, int depth=0 ) const;
			double get_sum_range( int index, int start=0, int stop=-1 );
			int get_count( int index ) const;
			int get_count() const;
			
			int get_dimension() const;
			
			bool next();
			void reset();
			
		private:
			looplist_t p_looplist;
			int p_count;
	};
}


#endif //util_h_flag

