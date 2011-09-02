//==========================================================================
// Title: util.cpp
//
// Author: Til Florian Guenzler & Christian G. Schroer
//
// description: adaptation to endian of machine
//==========================================================================

#include "util.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm> // for std::reverse
#include <cmath>
#include <ctime> // for std::time()
#include <cctype> // for std::tolower/tohigher
#ifndef WIN32
	#include <dirent.h>		// for reading directory entries with opendir()
	#include <sys/stat.h>	// for creating directories with mkdir()
	#include <termios.h>	// for set_keypress() and reset_keypress()
//	#include <unistd.h>	// for usleep()
#else
	// Set option '/clr' in order to make use of mscorlib.
	// Go to 'project properties:general:Use Managed Extensions'
	// You must then disable the Languague Extensions '/Za' option.
	#include <windows.h>
	#using <mscorlib.dll>	
	using namespace System::IO;
	using namespace System;
	#undef CreateDirectory // by JP on 21-11-2007

	// By JF:
	// Windows.h must be included, BEFORE using mscorlib, otherwise, strange errors will appear.
	// This way, Sleep() and matlibutil::usleep()can be used right away...
#endif

//***commented out the following 2 lines for the use in the cross correlator project, matlib needs this
//#include "constants.h"
//#include "commandline.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setprecision;
using std::setw;
using std::fstream;
using std::ifstream;
using std::getline;
using std::string;
using std::swap;
using std::vector;

#include <sstream>
using std::ostringstream;

//***commented out the following15 lines for the use in the cross correlator project, matlib needs this
//using ns_const::pi;
#define pi M_PI

namespace matlibutil
{
	const int MAX_COLS = 1000;
	const int MAX_ROWS = 1000;

	#ifndef WIN32
		const string TOKEN_START   			= "\033[";
		const string TOKEN_STOP				= "m";
		const string TOKEN_CLEAR_SCREEN_END	= "0J";
		const string TOKEN_CLEAR_SCREEN_BEG	= "1J";
		const string TOKEN_CLEAR_SCREEN_ALL	= "2J";
		const string TOKEN_CLEAR_LINE_END	= "0K";
		const string TOKEN_CLEAR_LINE_BEG	= "1K";
		const string TOKEN_CLEAR_LINE_ALL	= "2K";
		const string TOKEN_CURSOR_UP		= "A";
		const string TOKEN_CURSOR_DOWN		= "B";
		const string TOKEN_CURSOR_LEFT		= "D";
		const string TOKEN_CURSOR_RIGHT 	= "C";
		const string TOKEN_CURSOR_TO_COL	= "G";
		const string TOKEN_CURSOR_TO_POS	= "H";
		const string TOKEN_HIDE_CURSOR		= "?25l";
		const string TOKEN_SHOW_CURSOR		= "?25h";
		const string TEXT_RESET				= "0";
		const string BOLD_ON				= "1"; 
		const string BOLD_OFF				= "22";
		const string FAINT_ON				= "2"; 
		const string FAINT_OFF				= "22";
		const string ITALIC_ON				= "3"; 
		const string ITALIC_OFF 			= "23";
		const string UNDERLINE_1_ON  		= "4"; 
		const string UNDERLINE_2_ON  		= "21";
		const string UNDERLINE_OFF			= "24";
		const string BLINK_SLOW_ON  		= "5"; 
		const string BLINK_FAST_ON  		= "6"; 
		const string BLINK_OFF				= "25";
		const string INVERT_ON 				= "7"; 
		const string INVERT_OFF				= "27";
		const string HIDE_ON				= "8"; 
		const string HIDE_OFF				= "28";
		const string STRIKE_THRU_ON 		= "9"; 
		const string STRIKE_THRU_OFF		= "29";
	#else
		const string TOKEN_START			= "";
		const string TOKEN_CLEAR_LINE_END	= "";
		const string TOKEN_CLEAR_LINE_BEG	= "";
		const string TOKEN_CLEAR_LINE_ALL	= "";
		const string TOKEN_CLEAR_SCREEN_END	= "";
		const string TOKEN_CLEAR_SCREEN_BEG	= "";
		const string TOKEN_CLEAR_SCREEN_ALL	= "";
		const string TOKEN_CURSOR_UP		= "";
		const string TOKEN_CURSOR_DOWN		= "";
		const string TOKEN_CURSOR_LEFT		= "";
		const string TOKEN_CURSOR_RIGHT 	= "";
		const string TOKEN_CURSOR_TO_COL	= "";
		const string TOKEN_CURSOR_TO_POS	= "";
		const string TOKEN_HIDE_CURSOR		= "";
		const string TOKEN_SHOW_CURSOR		= "";
		const string TOKEN_STOP				= "";
		const string TEXT_RESET				= "";
		const string BOLD_ON			 	= "";
		const string BOLD_OFF				= "";
		const string FAINT_ON			 	= "";
		const string FAINT_OFF				= "";
		const string ITALIC_ON			 	= "";
		const string ITALIC_OFF 			= "";
		const string UNDERLINE_1_ON  	 	= "";
		const string UNDERLINE_2_ON  		= "";
		const string UNDERLINE_OFF			= "";
		const string BLINK_SLOW_ON  	 	= "";
		const string BLINK_FAST_ON  	 	= "";
		const string BLINK_OFF				= "";
		const string INVERT_ON 			 	= "";
		const string INVERT_OFF				= "";
		const string HIDE_ON			 	= "";
		const string HIDE_OFF				= "";
		const string STRIKE_THRU_ON 	 	= "";
		const string STRIKE_THRU_OFF		= "";
	#endif

	//--------------------------------------------------------------------------
	// Store the old terminal settings in this static variable.
	// It will be used by the functions 'set_keypress()' and 'reset_keypress()'.
	//--------------------------------------------------------------------------
#ifndef WIN32	//by JF: causes errors on vc7 -> use only on unix type systems
	static struct termios stored_settings;
#endif

	//==========================================================================
	//                                                               Simple math
	//==========================================================================

	//--------------------------------------------------------- placeValueSystem
	//
	// Calculates the representation of <value> with respect to place-value
	// system of given base. The representation is stored within <digits>
	// The first integer (digits[0]) will be the place with the lowest value,
	// the last integer (digits[digits.size()-1) will be the place with the
	// highest value, and it will be zero if and only if <value> is zero.
	// <digits> will be empty if <base> is smaller than two.
	// The sign of value will be ignored.
	//
	//--------------------------------------------------------------------------
	void placeValueSystem( double value, int base, vector<int> &digits )
	{
		vector<int>::size_type max_digits;
		
		digits.clear();
		value = floor( fabs( value ) );
		
		if ( value < 1 )
		{
			// value == 0
			digits.push_back( 0 );
			return;
		}
		
		// The maximal number of digits is given by the maximal size of
		// vector<int>.
		max_digits = digits.max_size();
		
		for ( vector<int>::size_type i=0; (i<max_digits) && (value >= pow(base,(int)i)); i++ )
		{
			digits.push_back( (int )floor( ( fmod( value, pow( base, (int)i+1 ) ) ) / pow( base, (int)i ) ) );
		}
	}

	//--------------------------------------------------------- placeValueSystem
	//
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
	//
	//--------------------------------------------------------------------------
	void placeValueSystem( double value, int base, int precision,
						   vector<int> &integral_digits,
						   vector<int> &fraction_digits )
	{
		double integral_part;
		double fraction_part;

		if ( precision > (int )fraction_digits.max_size() )
		{
			// This should never be true on most systems, because the maximal
			// size of vector<int> is 2^30-1 which is the highest value that
			// can be hold by integer variables.
			precision = (int )fraction_digits.max_size();
		}
																		// Example:
																		// value = -314.159265 | 41.996 | 0.0123
																		// base = 10
																		// precision = 2
																	
		integral_part = floor( fabs(value) );							// 314			| 41		| 0
		fraction_part = fabs(value) - integral_part;					// 0.159265 	| 0.996 	| 0.0123
		fraction_part *= pow( base, precision );						// 15.9265		| 99.6		| 1.23
		fraction_part = floor( fraction_part + 0.5 );					// 16			| 100.1 	| 1
		
		if ( fraction_part >= pow( base, precision ) )					// false		| true		| false
		{
			// The rounding affects the integral part.
			integral_part++;											//				| 42
			fraction_part = 0;											//				| 0
		}
		
		placeValueSystem( integral_part, base, integral_digits );		// <4, 1, 3>	| <2, 4>	| <0>
		
		if ( precision > 0 )
		{
			placeValueSystem( fraction_part, base, fraction_digits );	// <6, 1>		| <0>		| <1>
		}
		
		// Fill fractional digits with zeros until its size equals <precision>.
		for ( int i = (int )fraction_digits.size(); i<precision; i++ )
		{
			fraction_digits.push_back( 0 ); 							//	<6, 1>		| <0, 0>	| <1, 0>
		}

		// Reverse fractional_digits.
		if ( fraction_digits.size() > 0 )
		{
			std::reverse( fraction_digits.begin(), fraction_digits.end() );	// <1, 6>		| <0, 0>	| <0, 1>
		}
	}

	//----------------------------------------------------------------- getDigit
	//
	// Calculates the representation of <value> with respect to place-value
	// system of given base and returns the digit at given place.
	// Places of the integral part of the representation (the digits to the left
	// of the fractional point) are counted with positive numbers (starting with
	// zero for the less significant place) whereas place of the fractional part
	// of the representation are counted with negative number.
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
	//
	//--------------------------------------------------------------------------
	int getDigit( double value, int base, int place, int precision )
	{
		vector<int> integral_digits;
		vector<int> fraction_digits;
		int digit;
		
		if ( precision <= 0 )
		{		
			if ( place < 0 )
			{
				precision = -place + 1;
			}
			else
			{
				precision = 1;
			}
		}
		
		placeValueSystem( value, base, precision, integral_digits, fraction_digits );
		
		if ( place >= 0 )
		{
			if ( place < (int )integral_digits.size() )
			{
				digit = integral_digits[place];
			}
			else
			{
				// A (redundant) preceding zero.
				digit = 0;
			}
		}
		else
		{
			if ( (-place-1) < (int )fraction_digits.size() )
			{
				digit = fraction_digits[-place-1];
			}
			else
			{
				// A (redundant) successive zero.
				digit = 0;
			}
		}
		
		return digit;
	}

	//---------------------------------------------------------- newtonIteration
	//
	// Searches for a solution of f(x)=y using Newton iteration.
	// Starts at xstart and iterates as long as |f(x)-y|>|delta|.
	// g is the derivative of f.
	// The vector 'parameters' defines optional parameters.
	// Does only work if xstart is within the convergence interval.
	// (Should be improved with respect to this matter.)
	//
	//--------------------------------------------------------------------------
	double newtonIteration( newtonfunc_t f, newtonfunc_t g, const vector<double > &parameters,
							double y, double xstart, double delta )
	{
		double x;
		double fval;
		delta = fabs(delta);
		x = xstart;
		fval = f(x, parameters);
		while ( fabs( fval - y ) > delta )
		{
			x += ( y - fval ) / g(x, parameters);
			fval = f(x, parameters);
		}		
		return x;
	}

	//--------------------------------------------------------------------- pmod
	// Calculates the modulo of i with respect to n.
	// In contrast to the built-in operator '%', pmod will always have positive
	// results. While '-2%5' will result in '-2', 'pmod(-2,5)' will return '3'.
	//--------------------------------------------------------------------------
	int pmod( int i, int n )
	{
		return ((i%n)+n)%n;
	}


	//==========================================================================
	//                                                                Byte order
	//==========================================================================

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	int BigEndianByteOrder()
	{
    	short int word = 0x0001;
    	char *byte = (char *) &word;
    	return(byte[0] ? 0 : 1);
	}

	//--------------------------------------------------------- reverseByteOrder
	// Reverses the byte order of the integer numbers in a data field.
	// char *data: data field to be modified
	// size      : size of the datatype (in bytes)
	// length    : length of the data field
	// (by JP)
	//--------------------------------------------------------------------------
	void reverseByteOrder( char *data, int sizeDataType, int length )
	{
		for ( int i = 0; i < length; i++ )
		{
			for ( int j = 0; j < sizeDataType/2; j++ )
			{
				char temp;
				temp = data[i*sizeDataType+j];
				data[i*sizeDataType+j] = data[(i+1)*sizeDataType-j-1];
				data[(i+1)*sizeDataType-j-1] = temp;
			}
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(float *fval)
	{
    	float temp = *fval;
    	for (unsigned int i = 0; i < sizeof(float); i++)
    	{
        	((char *)(fval))[sizeof(float)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(double *dval)
	{
    	double temp = *dval;
    	for (unsigned int i = 0; i < sizeof(double); i++)
    	{
        	((char *)(dval))[sizeof(double)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(short *sval)
	{
    	short temp = *sval;
    	for (unsigned int i = 0; i < sizeof(short); i++)
    	{
        	((char *)(sval))[sizeof(short)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(int *ival)
	{
    	int temp = *ival;
    	for (unsigned int i = 0; i < sizeof(int); i++)
    	{
        	((char *)(ival))[sizeof(int)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(long *lval)
	{
    	long int temp = *lval;
    	for (unsigned int i = 0; i < sizeof(long); i++)
    	{
        	((char *)(lval))[sizeof(long)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(unsigned short *sval)
	{
    	unsigned short temp = *sval;
    	for (unsigned int i = 0; i < sizeof(unsigned short); i++)
    	{
        	((char *)(sval))[sizeof(unsigned short)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(unsigned int *ival)
	{
    	unsigned int temp = *ival;
    	for (unsigned int i = 0; i < sizeof(unsigned int); i++)
    	{
        	((char *)(ival))[sizeof(unsigned int)-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(unsigned long *lval)
	{
    	unsigned long temp = *lval;
    	for (unsigned int i = 0; i < sizeof(unsigned long ); i++)
    	{
        	((char *)(lval))[sizeof(unsigned long )-i-1] = ((char *)(&temp))[i];
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(float *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((float *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(double *data, unsigned long length)
	{
		for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((double *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(short *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((short *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(int *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((int *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(long *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((long *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(unsigned short *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((unsigned short *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(unsigned int *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((unsigned int *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void revByteOrder(unsigned long *data, unsigned long length)
	{
    	for (unsigned long i = 0; i < length; i++)
		{
			revByteOrder((unsigned long *)&(data[i]));
		}
	}

	//--------------------------------------------------------------------------
	// switches byteorder if big endian
	//--------------------------------------------------------------------------
	unsigned short convlitend( unsigned short &var )
	{
    	if ( BigEndianByteOrder() )
    	{
        	unsigned short temp=0;

        	temp = var << 8;
        	var >>= 8;
        	var += temp;
		}
       	return var;
	};


	//==========================================================================
	//                                                         vector operations
	//==========================================================================

	double scalarProd(vector<double> &v1, vector<double> &v2)
	{
		double value = 0.;
		if (v1.size() == v2.size())
		{
			for (unsigned int i = 0; i < v1.size(); i++)
				value += v1[i] * v2[i];
		}
		else
		{	
			cerr << "Error! Size of vectors does not match in scalarProd!" << endl;
		}
		return (value); 
	}

	//==========================================================================
	//                                                      random distributions
	//==========================================================================

	//--------------------------------------------------------------------------
	// equally distributed random numbers in range from low to high
	//--------------------------------------------------------------------------
	int equipart(int low, int high)
	{
    	return((rand()%(high-low+1))+low);
	}

	//--------------------------------------------------------------------------
	// equally distributed random numbers in range from low to high
	//--------------------------------------------------------------------------
	unsigned int equipart(unsigned int low, unsigned int high)
	{
    	return((rand()%((int )high-(int )low+1))+low);
	}

	//--------------------------------------------------------------------------
	// equally distributed random numbers in range from low to high
	//--------------------------------------------------------------------------
	double equipart(double low, double high)
	{
    	return(rand()*(high-low)/RAND_MAX+low);
	}

	//--------------------------------------------------------------------------
	// exponentially distributed random numbers
	//--------------------------------------------------------------------------
	int exprand(double E)
	{
    	int result = 0;
    	double y = equipart(0., 1.);
    	if (y > exp(E))
		result = 0;
    	else
		result = 1;
    	return(result);
	}

	//--------------------------------------------------------------------------
	// random points inside a sphere of radius 1 by von Neumann rejection
	//--------------------------------------------------------------------------
	int randsphere(double& x, double& y, double& z)
	{
		do {
			x = equipart(-1., 1.);
			y = equipart(-1., 1.);
			z = equipart(-1., 1.);
		} while (x*x+y*y+z*z > 1.);
		return(0);
	}

	/*
	//--------------------------------------------------------------------------
	// poisson noise according to J. Schnakenberg
	//--------------------------------------------------------------------------
	long int poisson(double average)
	{
    	if (average<0)
		average *= -1;
    	if (average < 100.)
    	{
		double c = exp(-average);
		long int n = 0;
		double Q = 1.;
		while (Q > c)
		{
	    	n++;
	    	Q *= (1.-(double )rand()/(double )RAND_MAX);
		}
		return(--n);
    	}
    	else
		return((long int )ceil(gauss(average)-0.5));
	}
	*/

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	double poisson(double average)
	{
		if( average == 0. )
			return 0.;
    	if( average < 0 )
			average *= -1;

    	if (average < 100.)
    	{
			double c = exp(-average);
			double n = 0.;
			double Q = 1.;
			while (Q > c)
			{
				n++;
				Q *= (1.-(double )rand()/(double )RAND_MAX);
			}
			return(--n);
    	}
    	else
		{
			double outGauss = gauss(average);
			if( outGauss < 0. )
				return 0.;
			else 
				return outGauss;
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void boxmueller(double& x1, double& x2)
	{
    	double y1, y2, s;
    	do
    	{
		y1 = 2*((double )rand()/(double )RAND_MAX) - 1.;
		y2 = 2*((double )rand()/(double )RAND_MAX) - 1.;
		s = y1*y1 + y2*y2;
    	} while (s >= 1.);

    	x1 = y1 * sqrt(-2.*log(s)/s);
    	x2 = y2 * sqrt(-2.*log(s)/s);
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void gauss(double& x1, double& x2, double average, double stddev)
	{
    	boxmueller(x1, x2);
    	x1 = average + x1*stddev;
    	x2 = average + x2*stddev;
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	double gauss(double average, double stddev)
	{
    	double x1, x2;
    	boxmueller(x1, x2);
    	return(average + x1*stddev);
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	double gausstrunc(double average, double stddev, double highlimit, double lowlimit)
	{
    	double x1, x2;
    	double value;
    	do {
			boxmueller(x1, x2);
			value = average + x1*stddev;
			if ((value > highlimit) || (value < lowlimit))
				value = average + x2*stddev;
    	} while ((value > highlimit) || (value < lowlimit));
    	return(value);
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	double gauss(double average)
	{
    	return(gauss(average, sqrt(average)));
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	int selectBin(vector<double> width)
	{
    	double sum = 0.;
    	for (unsigned int i = 0; i < width.size(); i++)
			sum += width[i];
    	if (sum == 0.)
			return(1);//Error: No selection available
    	double number = equipart(0., sum);
    	sum = 0.;
    	for (unsigned int i = 0; i < width.size(); i++)
    	{
			sum += width[i];
			if (sum > number)
				return (i);
    	}
    	return(0); // is never used
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	int bisect(vector<double> array) // 2D-array has to be monotonic!
	{
    	if (array.size() > 1)
    	{
		int low = 0;
		int up = (int)array.size();

		double p = equipart(array[0], array[array.size()-1]);
		//cout << "p: " << p << endl;

		while( (up-low) > 1 )
		{
	    	int mid = (low + up)/2;
	    	if(p >= array[mid]) low = mid;
	    	else up = mid;
	    	//cout << "low: " << low << ", up: " << up << endl;
		}
		return low;
    	} 
    	return 0;  
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void heapsort(vector<double> &array)
	{
    	unsigned int i, ir, j, l;
    	double tmp;

    	unsigned int n = (unsigned int)array.size();
    	if (n < 2) return;

    	l = n/2+1;
    	ir = n;

    	for(;;) 
    	{
		if( l>1 )
		{
	    	tmp = array[--l-1];
	    	cout << "l: " << l << endl;
		}
		else 
		{
	    	tmp = array[ir-1];
	    	array[ir-1] = array[0];
	    	if (--ir == 1)
	    	{
			array[0] = tmp;
			break;
	    	}
	    	cout << "ir: " << ir << endl;
		}
		i = l;
		j = l+l;

		while(j <= ir)
		{
	    	//cout << "array[" << j-1 << "]: " << array[j-1] << ", array[" << i-1 << "]: " << array[i-1] << endl;
	    	if(j < ir && (array[j-1] < array[j])) j++;
	    	if( tmp < array[j-1] ) 
	    	{
			array[i-1] = array[j-1];
			i = j;
			j *= 2; 
	    	} else j = ir+1;
		}
		array[i-1] = tmp;     
		//for (unsigned int j = 0; j < array.size(); j++)
		//{
		//  cout << array[j] << ' ';
		//}
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void indextable(vector<double> array, vector<unsigned int> &indexVec)
	{
    	indexVec.clear();
    	indexVec.resize(array.size());
    	for (unsigned int g = 0; g < indexVec.size(); g++)
			indexVec[g] = g;
		
    	unsigned int straight = 7;
    	unsigned int nstack = 50;
    	vector<unsigned int> istack;
    	istack.resize(nstack);
    	for (unsigned int h = 0; h < nstack; h++)
			istack.push_back(0);
    	unsigned int jstack = 0;
		
    	unsigned int i,j,k;
    	unsigned int ir = (unsigned int)array.size()-1;
    	unsigned int l = 0;
    	unsigned int itmp;
    	double vtmp;
		
    	for(;;)
    	{
			if (ir-l < straight) 
			{
				for( j = l+1; j <= ir; j++)
				{
					itmp = indexVec[j];
					vtmp = array[itmp];
					for ( i = j; i > l; i-- )
					{
						if (array[indexVec[i-1]] <= vtmp ) break;
						indexVec[i] = indexVec[i-1];
					}
					indexVec[i] = itmp;
				}
				
				if( jstack == 0 ) break;
				ir = istack[jstack--];
				l = istack[jstack--];
			}
			else 
			{
				k = (l+ir) >> 1;
				swap(indexVec[k], indexVec[l+1]);
				
				if (array[indexVec[l]] > array[indexVec[ir]])
					swap(indexVec[l],indexVec[ir]);
				if (array[indexVec[l+1]] > array[indexVec[ir]])
					swap(indexVec[l+1],indexVec[ir]);
				if (array[indexVec[l]] > array[indexVec[l+1]])
					swap(indexVec[l],indexVec[l+1]);
				
				i = l+1; 
				j = ir;
				itmp = indexVec[l+1];
				vtmp = array[itmp];
				
				for(;;)
				{
					do i++; while (array[indexVec[i]] < vtmp);
					do j--; while (array[indexVec[j]] > vtmp);
					
					if (j < i) break;
					swap(indexVec[i],indexVec[j]);
				}
				
				indexVec[l+1] = indexVec[j];
				indexVec[j] = itmp;
				jstack +=2;
				if (jstack > nstack) cerr << "NSTACK too small in indexing" << endl;
				
				if (ir-i+1 >= j-l)
				{
					istack[jstack] = ir;
					istack[jstack-1] = i;
					ir = j-1;
				}
				else 
				{
					istack[jstack] = j-1;
					istack[jstack-1] = l;
					l = i;
				}
			}
    	}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void init_rand()
	{
    	time_t now = std::time(0);
    	srand((unsigned int)now);
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	void isotrope(double& theta, double& phi)
	{
    	phi = equipart(0., 2.*3.1415926535);
    	theta = acos(equipart(-1., 1.));
	}


	//==========================================================================
	//                                                         string formatting
	//==========================================================================

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	string appendNumber(string base, unsigned int number, int num_digits)
	{
    	return(appendNumber(base, (int )number, num_digits));
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	string appendNumber(string base, long int number, int num_digits)
	{
    	return(appendNumber(base, (int )number, num_digits));
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	string appendNumber(string base, long unsigned  int number, int num_digits)
	{
    	return(appendNumber(base, (int )number, num_digits));
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	string appendNumber(string base, int number, int num_digits)
	{
 		string result("");
		result = base;
		ostringstream osst;

		osst << number;
		if ( num_digits > (int )osst.str().length() )
		{
			result.append( num_digits - osst.str().length(), '0' );
		}
		result += osst.str();
		return result;
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	string appendNumber(string base, double number, int length, int num_digits)
	{
    	ostringstream number_osst;
    	number_osst << setw(length);
    	number_osst << setprecision(num_digits);
    	number_osst << number;
    	string number_str( number_osst.str() );
	//	if ( string::size_type(number_digits) > number_str.length() )
	//		base.append( num_digits-number_str.length(), ' ' );
    	base += number_str;
    	return base;
	}

	//--------------------------------------------------------------------------
	// added on Nov,16,2005 by Jan Feldkamp
	//--------------------------------------------------------------------------
	string appendNumber(string left, string number, string right, string extension, int digits)
	{
		string ExtendedNumber = "";
		string preceedingZeros = "";
		int InitialDigits = (unsigned int)number.length();

		if (InitialDigits >= digits)
		{
			ExtendedNumber = number;
		}
		else
		{
			{
				for (int i = 0; i < (digits - InitialDigits); i++)
					preceedingZeros += "0";
				ExtendedNumber = preceedingZeros + number;
			}
		}

		return (left + ExtendedNumber + right + "." + extension);
	}

	//--------------------------------------------------------------------------
	// similar to 'generateFn()' but with c++ code instaed of c
	//--------------------------------------------------------------------------
	int appendNumber(string& base, string left, int number, string right, int number_digits)
	{
    	ostringstream number_osst;
    	number_osst << number;
    	string number_str( number_osst.str() );
    	base += left;
    	if ( string::size_type(number_digits) > number_str.length() )
		base.append( number_digits-number_str.length(), '0' );
    	base += number_str;
    	base += right;
    	return 0;
	}

	//----------------------------------------------------------------- toString
	//
	// Introduced for symmetry reasons.
	//
	//--------------------------------------------------------------------------
	string toString(string strval)
	{
    	return strval;
	}

	//----------------------------------------------------------------- toString
	//
	//--------------------------------------------------------------------------
	string toString(int number)
	{
    	ostringstream osst;
    	osst << number;
    	return osst.str();
	}

	//----------------------------------------------------------------- toString
	//
	//--------------------------------------------------------------------------
	string toString(int number, int digits)
	{
 		string result("");
		ostringstream osst;

		osst << number;
		if ( digits > (int) osst.str().length() )
		{
			result.append( digits - osst.str().length(), '0' );
		}
		result += osst.str();
		return result;
	}

	//----------------------------------------------------------------- toString
	//
	//--------------------------------------------------------------------------
	string toString(long int number)
	{
    	ostringstream osst;
    	osst << number;
    	return osst.str();
	}

	//----------------------------------------------------------------- toString
	//
	//--------------------------------------------------------------------------
	string toString(unsigned int number)
	{
    	ostringstream osst;
    	osst << number;
     	return osst.str();
	}

	//----------------------------------------------------------------- toString
	//
	//--------------------------------------------------------------------------
	string toString(unsigned long int number)
	{
    	ostringstream osst;
    	osst << number;
    	return osst.str();
	}

	//----------------------------------------------------------------- toString
	//
	//--------------------------------------------------------------------------
	string toString(double number, int prec)
	{
    	ostringstream osst;
    	osst << setprecision(prec) << number;
    	return osst.str();
	}


	//-------------------------------------------------------------------- toInt
	//
	//--------------------------------------------------------------------------
	int toInt(string number)
	{
		return (int)atol(number.c_str()); // calls cstdlib's atol (ascii to long int)
	}

	//----------------------------------------------------------------- toDouble
	//
	//--------------------------------------------------------------------------
	double toDouble(string number)
	{
		return (atof(number.c_str())); // calls cstdlib's atof (ascii to float)
	}

	//------------------------------------------------------------------ toLower
	//
	//--------------------------------------------------------------------------
	string toLower( string input )
	{
		string output;
		for ( string::iterator i=input.begin(); i!=input.end(); i++ )
		{
			output.append( 1, (char )std::tolower(*i) );
		}
		return output;
	}
	
	//------------------------------------------------------------------ toUpper
	//
	//--------------------------------------------------------------------------
	string toUpper( string input )
	{
		string output;
		for ( string::iterator i=input.begin(); i!=input.end(); i++ )
		{
			output.append( 1, (char)std::toupper(*i) );
		}
		return output;
	}

	//-------------------------------------------------------------------- equal
	//
	// Compares the strings A and B case insensitively.
	// Returns true, if A and B are equal apart from lower/upper case.
	//
	//--------------------------------------------------------------------------
	bool equal( string A, string B )
	{
		return toUpper(A) == toUpper(B);
	}

	//------------------------------------------------------------- spacepadding
	//
	// Returns a string, so that 'text+string' has the length given by width.
	// Returns the empty string, if text is longer than width.
	// As an alternative, the filling character can be changed to any ascii
	// character. This can be done with <fill>.
	//
	//--------------------------------------------------------------------------
	string spacepadding( string text, int width, char fill )
	{
		if ( width > (int )text.length() )
		{
			return string( width - (int )text.length(), fill );
		}
		else
		{
			return "";
		}
	}

	//--------------------------------------------------------------------------
	void set_keypress()
	{
#ifndef WIN32
    	struct termios new_settings;
    	tcgetattr( 0, &stored_settings );
    	new_settings = stored_settings;
    	new_settings.c_lflag &= ~ICANON;
    	new_settings.c_lflag &= ~ECHO;	// Turn off local echo
    	new_settings.c_cc[VTIME] = 0;
    	tcgetattr( 0, &stored_settings );
    	new_settings.c_cc[VMIN] = 1;
    	tcsetattr( 0, TCSANOW, &new_settings );
#endif
	}

	//--------------------------------------------------------------------------
	void reset_keypress()
	{
#ifndef WIN32
    	tcsetattr( 0, TCSANOW, &stored_settings );
#endif
	}

	//--------------------------------------------------------------------------
	int get_cursorPosX()
	{
		int row;
		int col;
		
		get_cursorPos( row, col );
		return col;
	}
	
	//--------------------------------------------------------------------------
	int get_cursorPosY()
	{
		int row;
		int col;
		
		get_cursorPos( row, col );
		return row;
	}

	//--------------------------------------------------------------------------
	void get_cursorPos( int &row, int &col )
	{
		row = 0;
		col = 0;
		
		set_keypress();
	    cout << "\033[6n";
    	scanf( "\033[%d;%dR", &row, &col );
		reset_keypress();
	}

	//--------------------------------------------------------------------------
	string clearLine()
	{
		return TOKEN_START + TOKEN_CLEAR_LINE_ALL;
	}

	//--------------------------------------------------------------------------
	string clearLineToBegin()
	{
		return TOKEN_START + TOKEN_CLEAR_LINE_BEG;
	}

	//--------------------------------------------------------------------------
	string clearLineToEnd()
	{
		return TOKEN_START + TOKEN_CLEAR_LINE_END;
	}

	//--------------------------------------------------------------------------
	string clearScreen()
	{
		return TOKEN_START + TOKEN_CLEAR_SCREEN_ALL;
	}

	//--------------------------------------------------------------------------
	string clearScreenToBegin()
	{
		return TOKEN_START + TOKEN_CLEAR_SCREEN_BEG;
	}

	//--------------------------------------------------------------------------
	string clearScreenToEnd()
	{
		return TOKEN_START + TOKEN_CLEAR_SCREEN_END;
	}

	//--------------------------------------------------------------------------
	string cursorUp( int num_lines )
	{
		if ( num_lines < 0 )
		{
			return cursorDown( -num_lines );
		}
		else
		{
			return TOKEN_START + toString( num_lines ) + TOKEN_CURSOR_UP;
		}
	}

	//--------------------------------------------------------------------------
	string cursorDown( int num_lines )
	{
		if ( num_lines < 0 )
		{
			return cursorUp( -num_lines );
		}
		else
		{
			return TOKEN_START + toString( num_lines ) + TOKEN_CURSOR_DOWN;
		}
	}

	//--------------------------------------------------------------------------
	string cursorLeft( int num_cols )
	{
		if ( num_cols < 0 )
		{
			return cursorRight( -num_cols );
		}
		else
		{
			return TOKEN_START + toString( num_cols ) + TOKEN_CURSOR_LEFT;
		}
	}

	//--------------------------------------------------------------------------
	string cursorRight( int num_cols )
	{
		if ( num_cols < 0 )
		{
			return cursorDown( -num_cols );
		}
		else
		{
			return TOKEN_START + toString( num_cols ) + TOKEN_CURSOR_RIGHT;
		}
	}

	//--------------------------------------------------------------------------
	string cursorToPosition( int row, int col )
	{
		return TOKEN_START + toString( row ) + ";" + toString( col ) + TOKEN_CURSOR_TO_POS;
	}

	//--------------------------------------------------------------------------
	string cursorToCol( int col )
	{
		return TOKEN_START + toString( col ) + TOKEN_CURSOR_TO_COL;
	}
	
	//--------------------------------------------------------------------------
	// Attention! This has not been tested.
	//
	string cursorToRow( int row )
	{
		int col;
		
		col = get_cursorPosX();
		return cursorToPosition( row, col );
	}
	
	//--------------------------------------------------------------------------
	string cursorLineBegin()
	{
		return TOKEN_START + toString( 1 ) + TOKEN_CURSOR_TO_COL;
	}

	//--------------------------------------------------------------------------
	string cursorLineEnd()
	{
		return TOKEN_START + toString( MAX_COLS ) + TOKEN_CURSOR_TO_COL;
	}

	//--------------------------------------------------------------------------
	string hideCursor()
	{
		return TOKEN_START + TOKEN_HIDE_CURSOR;
	}
	
	//--------------------------------------------------------------------------
	string showCursor()
	{
		return TOKEN_START + TOKEN_SHOW_CURSOR;
	}
	
	//--------------------------------------------------------------------------
	string textreset()
	{
		return ANSIescape( TEXT_RESET );
	}

	//--------------------------------------------------------------------------
	string boldon()
	{
		return ANSIescape( BOLD_ON );
	}

	//--------------------------------------------------------------------------
	string boldoff()
	{
		return ANSIescape( BOLD_OFF );
	}

	//--------------------------------------------------------------------------
	string fainton()
	{
		return ANSIescape( FAINT_ON );
	}

	//--------------------------------------------------------------------------
	string faintoff()
	{
		return ANSIescape( FAINT_OFF );
	}

	//--------------------------------------------------------------------------
	string italicon()
	{
		return ANSIescape( ITALIC_ON );
	}

	//--------------------------------------------------------------------------
	string italicoff()
	{
		return ANSIescape( ITALIC_OFF );
	}

	//--------------------------------------------------------------------------
	string underline1on()
	{
		return ANSIescape( UNDERLINE_1_ON );
	}

	//--------------------------------------------------------------------------
	string underline2on()
	{
		return ANSIescape( UNDERLINE_2_ON );
	}

	//--------------------------------------------------------------------------
	string underlineoff()
	{
		return ANSIescape( UNDERLINE_OFF );
	}

	//--------------------------------------------------------------------------
	string blinkslowon()
	{
		return ANSIescape( BLINK_SLOW_ON );
	}

	//--------------------------------------------------------------------------
	string blinkfaston()
	{
		return ANSIescape( BLINK_FAST_ON );
	}

	//--------------------------------------------------------------------------
	string blinkoff()
	{
		return ANSIescape( BLINK_OFF );
	}

	//--------------------------------------------------------------------------
	string inverton()
	{
		return ANSIescape( INVERT_ON );
	}

	//--------------------------------------------------------------------------
	string invertoff()
	{
		return ANSIescape( INVERT_OFF );
	}

	//--------------------------------------------------------------------------
	string hideon()
	{
		return ANSIescape( HIDE_ON );
	}

	//--------------------------------------------------------------------------
	string hideoff()
	{
		return ANSIescape( HIDE_OFF );
	}

	//--------------------------------------------------------------------------
	string strikethruon()
	{
		return ANSIescape( STRIKE_THRU_ON );
	}

	//--------------------------------------------------------------------------
	string strikethruoff()
	{
		return ANSIescape( STRIKE_THRU_OFF );
	}

	//--------------------------------------------------------------------------
	string ANSIescape( string code )
	{
		#ifndef WIN32
			return string( TOKEN_START ) + code + string( TOKEN_STOP );
		#else
			// Windows' shell does not support this kind of manipulation.
			return string( "" );
		#endif
	}

	//==========================================================================
	//                                                                      math
	//==========================================================================

	void complexSqrt(double &re, double &im)
	{
    	double r = sqrt(re*re+im*im);
    	double angle = acos(re/r);
    	if (asin(im/r) < 0.)
		angle = 2.*pi - angle;
    	re = sqrt(r) * cos(angle/2.);
    	im = sqrt(r) * sin(angle/2.);
	}

	// numrec routine
	double erfcc(double x)
	{
    	double t,z,ans;

    	z=fabs(x);
    	t=1.0/(1.0+0.5*z);
    	ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
					t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
					t*(-0.82215223+t*0.17087277)))))))));
    	if ( x >= 0 )
		{
			return ans;
		}
    	else 
		{
			return 2.0-ans;
		}
	}

	double erfkn(double x)
	{
    	return erfcc(-x)/2.;
	}

	double fermiFct(double x, double prec) 
	{
		double value = 0;
		if( prec < 1e-16 )
		{
			if( x < 0 )
				value = 1;
			else 
				value = 0;
		}
		else
		{
			double s = prec/2.197225; // =prec*2*ln3
			value = 1/(1+exp(x/s));
		}
		return value;
	}

	//==========================================================================
	//                                      Platform independant file management
	//==========================================================================
	
	//------------------------------------------------------------ get_separator
	//
	// Returns platform dependant path tree separator ( '/' or '\' ).
	//
	//--------------------------------------------------------------------------
	string get_separator()
	{
		#ifndef WIN32
			return "/";
		#else
			return "\\";
		#endif
	}

	//----------------------------------------------------- slashesToBackslashes
	//
	// Needed to convert path names from linux to windows convention.
	//
	//--------------------------------------------------------------------------
	string slashesToBackslashes( string path )
	{
		string output;

		for ( string::iterator i = path.begin(); i != path.end(); i++ )
		{
			if ( *i == '/' )
			{
				output += "\\";
			}
			else
			{
				output += ( *i );
			}
		}		
		return output;
	}

	//----------------------------------------------------- backslashesToSlashes
	//
	// Needed to convert path names from windows to linux convention.
	//
	//--------------------------------------------------------------------------
	string backslashesToSlashes( string path )
	{
		string output;

		for ( string::iterator i = path.begin(); i != path.end(); i++ )
		{
			if ( *i == '\\' )
			{
				output += "/";
			}
			else
			{
				output += ( *i );
			}
		}		
		return output;
	}

	//--------------------------------------------------------------- properPath
	//
	// Unix/Apple: convert backslashes to slashes.
	// MS Windows: convert slashes to backslashes.
	//
	//--------------------------------------------------------------------------
	string properPath( string path )
	{
		string output;

		#ifndef WIN32
			output = backslashesToSlashes( path );
		#else
			output = slashesToBackslashes( path );
		#endif

		return output;
	}

	//---------------------------------------------------- returnDirectoryPartOf
	//
	// Return the directory part of a path.
	// Works with slashes or backslashes.
	// Example: "this/is/it" -> "this/is/"
	// Example: "this\is\it" -> "this\is\"
	//
	//--------------------------------------------------------------------------
	string returnDirectoryPartOf( string path, bool verbose )
	{
		string output;
		string::size_type pos_sep;
		
		if (verbose)
		{
			cout << "\nreturnDirectoryPartOf(" << path << ")" << endl;
		}

		pos_sep = path.find_last_of( "/\\" );
/*
		pos_sep_slash = path.rfind( "/" );
		pos_sep_backslash = path.rfind( "\\" );


		if ( (pos_sep_slash == string::npos) && (pos_sep_backslash == string::npos) )
		{
			output = "";
		}
		else if ( (pos_sep_slash != string::npos) && (pos_sep_backslash != string::npos) )
		{
			string::size_type max_pos;
			
			if ( pos_sep_slash > pos_sep_backslash )
			{
				output = path.substr( 0, pos_sep_slash+1 );
			}
			else
			{
				output = path.substr( 0, pos_sep_backslash+1 );
			}
		}
		else if ( pos_sep_slash != string::npos )
		{
			output = path.substr( 0, pos_sep_slash+1 );		
		}
*/
		if (pos_sep == string::npos)
		{ 
			output = "";
		}
		else
		{
			output = path.substr( 0, pos_sep+1 );
		}
		if (verbose)
		{
			cout << "directory part: '" << output << "'" << endl; 
		}
		
		return output;
	}

	//----------------------------------------------------- returnFilenamePartOf
	//
	// Return the filename part of a path.
	// Example: this/is/it -> it
	// Example: this\is\it -> it
	//
	//--------------------------------------------------------------------------
	string returnFilenamePartOf( string path, bool verbose )
	{
		string output;
		string::size_type pos_sep;

		if (verbose)
		{
			cout << "\nreturnFilenamePartOf(" << path << ")" << endl;
		}

		pos_sep = path.find_last_of( "/\\" );

		if ( pos_sep != string::npos )
		{
			output = path.substr( pos_sep+1 );
		}
		else
		{
			output = path;
		}

		if (verbose)
		{
			cout << "Returning file name part: '" << output << "'" << endl; 
		}
		
		return output;
	}

	//---------------------------------------------------------- resolvePathTree
	//
	// Extracts the parts from a path tree and returns them within a vector of
	// strings. An absolute path results in a vector with the empty string as
	// the first element. Consecutive separators are treated as just one separator.
	// Example: "aaa/bbb/ccc/ddd" ---> ( "aaa", "bbb", "ccc", "ddd" )
	// Example: "/aaa/bbb/ccc"    ---> ( "", "aaa", "bbb", "ccc" )
	// Example: "aaa///bbb/ccc/"  ---> ( "aaa", "bbb", "ccc" )
	//
	//--------------------------------------------------------------------------
	int resolvePathTree( string path, vector<string> &parts )
	{
		return resolvePathTree( path, parts, *(get_separator().begin()) );
	}
	int resolvePathTree( string path, vector<string> &parts, char separator )
	{
		string part("");
		bool again;

		parts.clear();
		again = false;
		for ( string::iterator i=path.begin(); i!=path.end(); i++ )
		{
			if ( ( *i == separator ) && !again )
			{
				parts.push_back( part );
				part = "";
				again = true;
			}
			else if ( *i != separator )
			{
				part.push_back( *i );
				again = false;
			}
		}

		if ( part != "" )
		{
			parts.push_back( part );
		}

		return 0;
	}

	//----------------------------------------------------------- createPathTree
	//
	// Creates the full path tree, even if more than the last part of the path
	// does not exist. Fails, if any of the parts of the path tree already
	// exists as a non-directory file type.
	// Example: 'aaa/bbb' + createPathTree( "aaa/bbb/ccc/ddd" ) --> 'aaa/bbb/ccc/ddd'
	// Example: 'aaa/bbb' + createPathTree( "aaa" ) --> 'aaa/bbb'
	//
	//--------------------------------------------------------------------------
	int createPathTree( string pathtree, bool verbose )
	{
		vector<string> parts;
		string path;
		int i_start = 0;

		if ( verbose>3 )
		{
			cout << "pathtree: '" << pathtree << "'" << endl;
		}

		if ( pathtree == "" )
		{
			if ( verbose>2 )
			{
				cout << "Given pathtree is empty. Nothing done." << endl;
			}
			return 1;
		}

		#ifndef WIN32
			resolvePathTree( pathtree, parts );

			if ( verbose>2 )
			{
				cout << "parts.size()=" << parts.size() << endl;
			}

			if ( parts.size() == 0 )
			{
				if ( verbose>2 )
				{
					cout << "The vector 'parts' is empty. Nothing done." << endl;
				}
				return 2;
			}

			if ( parts[0] == "" )
			{
				if ( verbose>2 )
				{
					cout << "Absolute path." << endl;
				}
				path = "/";
				i_start = 1;
			}
			else
			{
				if ( verbose>2 )
				{
					cout << "relative path." << endl;
				}
				path = "";
				i_start = 0;
			}

			for ( int i=i_start; i<(int )parts.size(); i++ )
			{
				int retval;

				path += parts[i] + get_separator();
				if ( verbose>2 )
				{
 					cout << "part number " << i << ": '" << path << "'" << std::flush;
				}
				
				if ( !directoryDoesExist( path ) )
				{
					if ( verbose>1 )
					{
						cout << " <-- Directory does not exist. Create it now." << endl;
					}
					
					retval = mkdir( path.c_str(), 511 );
					if ( retval )
					{
						if ( verbose>0 )
						{
							cout << "ERROR. Could not create directory." << endl;
						}
						return 3;
					}
				}
			}
		#else
			if ( !directoryDoesExist( path ) )
			{	
				System::IO::Directory::CreateDirectory( path.c_str() );
			}
		#endif
		
		return 0;
	}


	//------------------------------------------------------------ fileDoesExist
	//
	// Returns true if 'fullname' does exist as a regular file.
	// Returns false if 'fullname' does not exist or if 'fullname'
	// is not a regular file (e.g. a directory or a fifo).
	//
	//--------------------------------------------------------------------------
	bool fileDoesExist( string fullname )
	{
		bool existing=0;

		#ifndef WIN32
			int retval;
			struct stat sb;

			retval = lstat( fullname.c_str(), &sb );

			if ( retval == -1 )
			{
				existing = false;
			}
			else
			{
				if ( S_ISREG( sb.st_mode ) )
				{
					existing = true;
				}
				else
				{
					// Exists but is not a regular file.
					existing = false;
				}
			}
		#else
			if ( System::IO::File::Exists( fullname.c_str() ) )
			{
				existing = true;
			}
		#endif

		return existing;
	}

	//------------------------------------------------------------ fifoDoesExist
	//
	// Returns true if 'fullname' does exist as a fifo file.
	// Returns false if 'fullname' does not exist or if 'fullname'
	// is not a fifo file (e.g. a directory or a regular file).
	//
	//--------------------------------------------------------------------------
	bool fifoDoesExist( string fullname )
	{
		bool existing=0;

		#ifndef WIN32
			int retval;
			struct stat sb;

			retval = lstat( fullname.c_str(), &sb );

			if ( retval == -1 )
			{
				existing = false;
			}
			else
			{
				if ( S_ISFIFO( sb.st_mode ) )
				{
					existing = true;
				}
				else
				{
					// Exists but is not a regular file.
					existing = false;
				}
			}
		#else
			// Alternative could be to always return false on WIN32.
			if ( System::IO::File::Exists( fullname.c_str() ) )
			{
				existing = true;
			}
		#endif

		return existing;
	}

	//------------------------------------------------------- directoryDoesExist
	//
	//	Returns true if 'path' does exist as a directory.
	//	Returns false if 'fullname' does not exist or if 'fullname'.
	//	is not a directory.
	//
	//--------------------------------------------------------------------------
	bool directoryDoesExist( string path, bool verbose )
	{
		bool existing=0;
		//bool verbose=1;

		if ( verbose )
		{
			cout << "Directory to check: '" << path << "'" << endl;
		}
		
		#ifndef WIN32								
			int retval;
			struct stat sb;

			retval = lstat( path.c_str(), &sb );
			
			if ( verbose )
			{
				cout << "lstat returned: " << retval << endl;
			}

			if ( retval == -1 )
			{
				// lstat failed.
				if ( verbose )
				{
					cout << "WARNING. lstat failed. Assume not existing." << endl;
				}
				existing = false;
			}
			else
			{
				if ( S_ISDIR( sb.st_mode ) )
				{
					existing = true;
				}
				else
				{
					// Exists but is not a directory.
					existing = false;
				}
			}
		#else
			existing = System::IO::Directory::Exists( path.c_str() );
		#endif

		if ( verbose )
		{
			cout << "Given directory does " << (existing?"":"not ") << "exist." << endl;
		}

		return existing;
	}

	//------------------------------------------------------- checkFileExistance
	//
	// Checks, if the file exists in the current directory.
	// Returns true if file exists and can be opened.
	// In contrast to 'fileDoesExist()', this functions fails,
	// if the given file exists but cannot be opened because
	// of missing permissions or something similar.
	//
	//--------------------------------------------------------------------------
	bool checkFileExistance( string FileName )
	{
		bool fileExists;
		fstream test;
				
		fileExists = false;
		test.open( FileName.c_str() );
		test.close();
		if ( test.fail() )
		{
			fileExists = false;
		}
		else
		{
			fileExists = true;
		}
		test.clear();
		return fileExists;
	}

	//--------------------------------------------------------------- countLines
	//
	// Returns the number of lines in a file.
	//
	//--------------------------------------------------------------------------
	long int countLines(string fn)
	{
    	long int linecount;
    	long int countcr;
    	long int countlf;
		ifstream ein;
    	
		linecount = -1;
    	countcr = 0;
    	countlf = 0;

    	ein.open( fn.c_str() );
    	if ( ein.good() )
    	{
			char temp;
			
			while ( !ein.eof() )
			{
	    		ein.get( temp );
	    		if ( temp == 13 )
				{
					countcr++;
				}
	    		if ( temp == 10 )
				{
					countlf++;
				}
			}
			
			if ( temp == 13 )
			{
				countcr--;
			}
			
			if ( temp == 10 )
			{
				countlf--;
			}
			
			cout << countcr << " " << countlf << endl;
			
			if ( countcr > countlf )
			{
	    		linecount = countcr;
			}
			else
			{
	    		linecount = countlf;
			}
    	}
    	return linecount;
	}

	//--------------------------------------------------------------- countLines
	//
	// Returns the number of lines in a file.
	//
	//--------------------------------------------------------------------------
	long int countLines(string fn, string skipchar)
	{
    	long int linecount;
    	long int countcr;
    	long int countlf;
    	ifstream ein;

    	linecount = -1;
    	countcr = 0;
    	countlf = 0;
    	ein.open( fn.c_str() );

    	if ( ein.good() )
    	{
			char temp;
			int newline = 1;
			int skipflag = 0;
			
			while ( !ein.eof() )
			{
	    		ein.get( temp );
				
	    		if ( newline )
	    		{
					skipflag = 0;
					newline = 0;
					
					skipflag = ( skipchar.find( temp ) != string::npos );
					/*
					for (unsigned int i = 0; i < skipchar.size(); i++)
					{
		    			if (temp == skipchar[i])
		    			{
							skipflag = 1;
		    			}
					}
					*/
	    		}

	    		if (temp == 13)
	    		{
					if (!skipflag)
					{
		    			countcr++;
					}
					newline = 1;
	    		}

	    		if (temp == 10)
	    		{
					if (!skipflag)
					{
		    			countlf++;
					}
					newline = 1; 
	    		}
			}
			
			if (temp == 13)
			{
				countcr--;
			}
			
			if (temp == 10)
			{
				countlf--;
			}
			
			cout << countcr << " " << countlf << endl;
			
			if (countcr > countlf)
			{
	    		linecount = countcr;
			}
			else
			{
	    		linecount = countlf;
			}
    	}
    	return linecount;
	}

	//==========================================================================
	//                                                               portability
	//==========================================================================

	//------------------------------------------------------------------- usleep
	int usleep( useconds_t time_usec ) 
	{
		#ifdef WIN32
			// define usleep() on windows machines so that it can be used
			// consistently for WIN32 and linux by specifying the namespace
			// matlibutil whenever calling usleep;
			// See note above, windows.h must be include, BEFORE using mscorlib.
			Sleep( (DWORD )ceil( (double )time_usec/1000.0 ) );
		#else
			// Unix and Mac OS X have the 'usleep()' routine built-in.
			// Important note.
			// Do not omit the double colon. Infinite recursive calling of
			// functions may be annoying...
			::usleep(time_usec);
		#endif

		return 0;
	}

	//==========================================================================
	//                                                                statistics
	//==========================================================================

	//-------------------------------------------------------------------- binio
	binio::binio()
	{
	}
	binio::~binio()
	{
	}
	int binio::write(FILE *fp, char *buffer, int size)
	{
		for (int i = 0; i < size; i++)
			putc(buffer[i], fp);
		return(0);
	}
	int binio::read(FILE *fp, char *buffer, int size)
	{	
		for (int i = 0; i < size; i++)
			buffer[i] = (char)getc(fp);
		return(0);
	}


	//---------------------------------------------------------------------- Min
	Min::Min()
	{
		reset();
	}
	Min::~Min()
	{
	}
	void Min::reset()
	{
		minimum = 1e100;
	}
	double Min::get()
	{
		return(minimum);
	}
	int Min::compare(double x)
	{
		if (minimum > x)
		{
			minimum = x;
			return 1;
		}
		else
		{
			return 0;
		}
	};


	//---------------------------------------------------------------------- Max
	Max::Max()
	{
		reset();
	}
	Max::~Max()
	{
	}
	void Max::reset()
	{
		maximum = -1e100;
	}
	double Max::get()
	{
		return maximum;
	}
	int Max::compare(double x)
	{
		if (maximum < x)
		{
			maximum = x;
			return 1;
		}
		else
		{
			return 0;
		}
	 }

	//-------------------------------------------------------------------- Range
	Range::Range()
	{
		reset();
	}
	Range::~Range()
	{
	}
	void Range::reset()
	{
		maximum = -1e100;
		minimum = 1e100;
	}
	double Range::get()
	{
		if (maximum >= minimum)
			return maximum - minimum;
		else 
			return 0.;
	}
	double Range::getMax()
	{
		return maximum;
	}
	double Range::getMin()
	{
		return minimum;
	}
	int Range::compare(double x)
	{
		int flag = 0;
		if (maximum < x)
		{
			maximum = x;
			flag = 1;
		}
		if (minimum > x)
		{
			minimum = x;
			flag = 1;
		}			
		return flag;
	}
		
	//--------------------------------------------------------------------- IMin
	IMin::IMin()
	{
		reset();
	}
	IMin::~IMin()
	{
	}
	void IMin::reset()
	{
		minimum = ((long int ) 1 << ((sizeof(long int) * 8) - 2));
	}
	long int IMin::get()
	{
		return minimum;
	}
	int IMin::compare(long int x)
	{
		if (minimum > x)
		{
			minimum = x;
			return 1;
		}
		else
		{
			return 0;
		}
	};


	//--------------------------------------------------------------------- IMax
	IMax::IMax()
	{
		reset();
	}
	IMax::~IMax()
	{
	}
	void IMax::reset()
	{
		maximum = ((long int ) 1 << ((sizeof(long int) * 8) - 1));
	}
	long int IMax::get()
	{
		return maximum;
	}
	int IMax::compare(long int x)
	{
		if (maximum < x)
		{
			maximum = x;
			return 1;
		}
		else
		{
			return 0;
		}
	};


	//-------------------------------------------------------------------- UIMax
	UIMax::UIMax()
	{
		reset();
	}
	UIMax::~UIMax()
	{
	}
	void UIMax::reset()
	{
		maximum = 0;
	}
	unsigned long int UIMax::get()
	{
		return(maximum);
	}
	int UIMax::compare(unsigned long int x)
	{
		if (maximum < x)
		{
			maximum = x; return(1);
		}
		else
		{
			return(0);
		}
	};


	//---------------------------------------------------------------------- sum
	sum::sum()
	{
		reset();
	}
	sum::~sum()
	{
	}
	void sum::reset()
	{
		count = 0;
		summe = 0.;
	}
	double sum::get()
	{
		return summe;
	}
	long int sum::getCount()
	{
		return count;
	}
	void sum::add(double x)
	{
		summe += x;
		count++;
	}

	//--------------------------------------------------------------------- mean
	mean::mean()
	{
    	reset();
	}
	mean::~mean()
	{
	}
	void mean::reset()
	{
    	count = 0; 
    	sum = 0.; 
	}
	double mean::get()
	{
    	if (count)
		{
			return sum/count;
		}
    	else
		{
			return 0;
		}
	}
	double mean::getSum()
	{
    	return sum;
	}
	long int mean::getCount()
	{
    	return count;
	}
	void mean::add(double x, unsigned int n)
	{
    	sum += x*n;
    	count += n;
	}


	//--------------------------------------------------------------------- histogram
	histogram::histogram(double lowval, double highval, int intervals)
	{
		// by JP
		if ( lowval >= highval )
		{
			cout << "WARNING. In histogram::histogram() - lowval must be smaller than highval."
				 << "Force highval to lowval+1." << endl;
			highval = lowval+1;
		}
		
		setHigh(highval);
		setLow(lowval);
		hist.resize(intervals);
		clear();
	}
	
	histogram::histogram(histogram *histoval)
	{
		double lowval;
		double highval;
		
		highval = histoval->getHigh();
		lowval = histoval->getLow();
		
		// by JP
		if ( lowval >= highval )
		{
			cout << "WARNING. In histogram::histogram() - lowval must be smaller than highval."
				 << "Force highval to lowval+1." << endl;
			highval = lowval+1;
		}
		
		setHigh(highval);
		setLow(lowval);
		hist.resize(histoval->getSize());
		clear();
		add(histoval);
	}
	
	void histogram::add(histogram *histoval)
	{
		if ((histoval->getHigh() == high) && (histoval->getLow() == low) && (histoval->getSize() == hist.size()))
		{
			for (unsigned int i = 0; i < hist.size(); i++)
			{
				hist[i] += histoval->get(i);
				count += histoval->get(i);
			}
		}
		else
		{
			cerr << "Error! Histograms do not share same range in add! No action performed!" << endl;
		}
	}
	
	// Add given value to the corresponding channel.
	void histogram::add(double val)
	{
		int channel;
		
		if ((val < high) && (val >=low))
		{
			channel = (unsigned int ) floor((val-low)/(high-low)*hist.size());
			hist[channel]++;
			count++;
		}
		else
		{
			cerr << "WARNING. In histogram::add() - " << val << " not a valid value in histogram" << endl;
		}
	}
	
	// Add to given channel
	void histogram::add(const unsigned int channel)
	{
		if ( channel < getSize() )
		{
			hist[channel]++;
			count++;
		}
		else
		{
			cerr << "WARNING. In histogram::add() - " << channel << " not a valid channel in histogram" << endl;
		}
	}
	
	void histogram::add_nocheck(const unsigned int channel)
	{
		hist[channel]++;
		count++;
	}
	
	void histogram::clear()
	{
		count = 0;
		for (unsigned int i = 0; i < hist.size(); i++)
		{
			hist[i] = 0;
		}
	}
	
	// Return number of counts within given channel.
	unsigned long int histogram::get(const unsigned int channel)
	{
		if ( channel < hist.size() )
		{
			return hist[channel];
		}
		else
		{	
			cerr << "WARNING. In histogram::get() - Histogram channel " << channel << " out of bounds." << endl;
			return 0;
		}
	}
	
	// Return number of counts within channel corresponding to given value.
	unsigned long int histogram::get(double val)
	{
		int channel = (int )floor( (val-low)/(high-low) );
		
		if ( ( channel >= 0 ) && ( channel < (int )hist.size() ) )
			return get( (unsigned int )channel );
		else
		{
			cerr << "WARNING. In histogram::get() - " << "Value " << val << " out of bounds." << endl;
			return(0);
		}
	}
	
	void histogram::get(vector<unsigned long int > &histval)
	{
		histval = hist;
	}

	void histogram::get(vector<double > &histval)
	{
		histval.resize(hist.size());
		for (unsigned int i = 0; i < hist.size(); i++)
		{
			histval[i] = hist[i];
		}
	}

	double histogram::getHigh()
	{
		return(high);
	}
	
	double histogram::getLow()
	{
		return(low);
	}
	
	unsigned long int histogram::getSize()
	{
		return( (unsigned long int)hist.size() );
	}
	
	double histogram::getMean(int center)
	{
		matlibutil::mean mittel;
		for (unsigned int i = 0; i < hist.size(); i++)
		{
			mittel.add((i+0.5*center)*(high-low)/hist.size()+low, (unsigned int ) get(i));
		}
		return(mittel.get());
	}
	
	double histogram::getStDev()
	{
		matlibutil::stDev stdev;
		for (unsigned int i = 0; i < hist.size(); i++)
		{
			stdev.add((i+0.5)*(high-low)/hist.size()+low, (unsigned int )get(i));
		}
		return(stdev.get());
	}

	double histogram::getErrorOfMean()
	{
		matlibutil::stDev stdev;
		for (unsigned int i = 0; i < hist.size(); i++)
		{
			stdev.add((i+0.5)*(high-low)/hist.size()+low, (unsigned int )get(i));
		}
		return(stdev.getErrorOfMean());
	}

	double histogram::getPercentile(double percentval)
	{
		if ((percentval >= 0.) && (percentval <=1.))
		{
			long int population = (long int )floor((double )count*percentval);
			if (population > 0)
			{
				matlibutil::sum summe;
				unsigned int i = 0; 
				while ((summe.get() < population) && (i < hist.size()))
				{
					summe.add(get(i));
					i++;
				}
				long int difference = (long int )floor(population - summe.get());
				if (i == hist.size())
					return (high);
				else
					return ((i+(double )difference/(double )get(i-1))*(high-low)/hist.size()+low);
			}
			else
				return(low);
		}
		else
		{
			cerr << "Error! Percentage out of range in getPercentile! Must be between 0 and 1!" << endl;
			return(0.);
		}
	}
	
	double histogram::getMedian()
	{
		return(getPercentile(0.5));
	}
	
	double histogram::getErrorOfMedian()
	{
		return(getErrorOfMedian(1000));
	}
	
	double histogram::getErrorOfMedian(unsigned int trials)
	{
		matlibutil::stDev stdev;
		vector<double > distribution;
		get(distribution);

		for (unsigned int i = 0; i < trials; i++)
		{
			histogram newdistribution(this);
			newdistribution.clear();
			for (unsigned int j = 0; j < getCount(); j++)
			{
				newdistribution.add_nocheck((int )selectBin(distribution));
			}
			stdev.add(newdistribution.getMedian());
		}
		return(stdev.get());
	}

	void histogram::setHigh(double highval)
	{
		high = highval;
	}
	
	void histogram::setLow(double lowval)
	{
		low = lowval;
	}
	
	unsigned long int histogram::getCount()
	{
		return(count);
	}
	
	// If consistent, this function should return the value stored in 'count'.
	unsigned long int histogram::sum()
	{	
		matlibutil::sum summe;
		for (unsigned int i = 0; i < hist.size(); i++)
			summe.add(hist[i]);
		return ((unsigned long int)summe.get());
	}


	//------------------------------------------------------------------ qMoment
	qMoment::qMoment()
	{
		reset(); q = 2.;
	}
	qMoment::~qMoment()
	{
	}
	void qMoment::setMoment(long int countval, double momentval)
	{
		count = countval;
		sum = momentval*countval;
	}
	void qMoment::addMoment(long int countval, double momentval)
	{
		count += countval;
		sum += momentval*countval;
	}
	int qMoment::setQ(double qval)
	{
		q = qval;
		return 0;
	}
	void qMoment::reset()
	{
		count = 0; sum = 0.;
	}
	double qMoment::get()
	{
		return pow(sum/count,1/q);
	}
	double qMoment::getQPow()
	{
		return sum/count;
	}
	long int qMoment::getCount()
	{
		return count;
	}
	void qMoment::add(double x)
	{
		sum += pow(x,q);
		count++;
	}

	//------------------------------------------------------------------ counter
	counter::counter()
	{
		reset();
	}
	counter::~counter()
	{
	}
	void counter::reset()
	{
		count = 0;
	}
	long int counter::getCount()
	{
		return count;
	}
	void counter::add()
	{
		count++;
	}


	//--------------------------------------------------------------------- Flag
	Flag::Flag(double comparatorval)
	{
		setComparator(comparatorval);
		reset();
	}
	Flag::~Flag()
	{
	}
	void Flag::reset()
	{
		set(0);
	}
	void Flag::setComparator(double comparatorval)
	{
		comparator = comparatorval;
	}
	double Flag::getComparator()
	{
		return comparator;
	}
	void Flag::set(int flagval)
	{
		flag = flagval;
	}
	int Flag::get()
	{
		return flag;
	}
	double Flag::compare(double x)
	{
		if (x == comparator)
		{
			flag = 1;
		}
	 	return x;
	 }


	//--------------------------------------------------------------------- Stat
	Stat::Stat() : sum(0.), quad(0.), count(0)
	{
	}
	Stat::~Stat()
	{
	}
	void Stat::set(Stat *statval)
	{
		quad = statval->getQuad();
		sum = statval->getSum(); 
		count = statval->getCount();
	}
	void Stat::add(Stat *statval)
	{
		quad += statval->getQuad();
		sum += statval->getSum(); 
		count += statval->getCount();
	}
	void Stat::reset()
	{
		count = 0;
		sum = 0;
		quad = 0;
	}
	double Stat::getStd()
	{
		return sqrt(quad/count-sum*sum/count/count);
	}
	double Stat::getStd2()
	{
		return quad/count-sum*sum/count/count;
	}
	double Stat::getMean()
	{
		return sum/count;
	}
	long int Stat::getCount()
	{
		return count;
	}
	void Stat::add(double x)
	{
		sum += x;
		count++;
		quad += x*x;
	}
	double Stat::getQuad()
	{
		return quad;
	}
	double Stat::getSum()
	{
		return sum;
	}

	//-------------------------------------------------------------------- statf
	statf::statf() : sum(0.), quad(0.), count(0)
	{
	}
	statf::~statf()
	{
	}
	void statf::set(statf *statval)
	{
		quad = statval->getQuad();
		sum = statval->getSum(); 
		count = statval->getCount(); }
	void statf::add(statf *statval)
	{
		quad += statval->getQuad();
		sum += statval->getSum(); 
		count += statval->getCount();
	}
	void statf::reset()
	{
		count = 0;
		sum = 0.;
		quad = 0.;
	}
	float statf::getStd()
	{
		return sqrt(quad/count-sum*sum/count/count);
	}
	float statf::getStd2()
	{
		return quad/count-sum*sum/count/count;
	}
	float statf::getMean()
	{
		return sum/count;
	}
	long int statf::getCount()
	{
		return count;
	}
	void statf::add(double x)
	{
		sum += (float)x;
		count++;
		quad += float(x*x);
	}
	float statf::getQuad()
	{
		return(quad);
	}
	float statf::getSum()
	{
		return(sum);
	}


	//----------------------------------------------------------------- twodstat
	twodstat::twodstat() : sum1(0.), sum2(0.), quad(0.), count(0)
	{
	}
	twodstat::~twodstat()
	{
	}
	void twodstat::set(twodstat *statval)
	{
		quad = statval->getQuad(); 
		sum1 = statval->getSum1();
		sum2 = statval->getSum2(); 
		count = statval->getCount();
	}
	void twodstat:: add(twodstat *statval)
	{
		quad += statval->getQuad(); 
		sum1 += statval->getSum1();
		sum2 += statval->getSum2(); 
		count += statval->getCount();
	}
	void twodstat::reset()
	{
		count = 0;
		sum1 = 0.;
		sum2 = 0.;
		quad = 0.;
	}
	double twodstat::getStd()
	{
		return sqrt(quad/count-(sum1*sum1+sum2*sum2)/count/count);
	}
	double twodstat::getStd2()
	{
		return quad/count-(sum1*sum1+sum2*sum2)/count/count;
	}
	double twodstat::getMeanX()
	{
		return sum1/count;
	}
	double twodstat::getMeanY()
	{
		return sum2/count;
	}
	long int twodstat::getCount()
	{
		return count;
	}
	void twodstat::add(double x, double y)
	{
		sum1 += x;
		sum2 += y;
		count++; 
		quad += x*x + y*y;
	}
	double twodstat::getSum1()
	{
		return sum1;
	}
	double twodstat::getSum2()
	{
		return sum2;
	}
	double twodstat::getQuad()
	{
		return quad;
	}


	//---------------------------------------------------------------- twodstatf
	twodstatf::twodstatf() : sum1(0.), sum2(0.), quad(0.), count(0)
	{
	}
	twodstatf::~twodstatf()
	{
	}
	void twodstatf::set(twodstatf *statval)
	{
		quad = statval->getQuad(); 
		sum1 = statval->getSum1();
		sum2 = statval->getSum2(); 
		count = statval->getCount();
	}
	void twodstatf::add(twodstatf *statval)
	{
		quad += statval->getQuad(); 
		sum1 += statval->getSum1();
		sum2 += statval->getSum2(); 
		count += statval->getCount();
	}
	void twodstatf::reset()
	{
		count = 0;
		sum1 = 0.;
		sum2 = 0.;
		quad = 0.;
	}
	float twodstatf::getStd()
	{
		return sqrt(quad/count-(sum1*sum1+sum2*sum2)/count/count);
	}
	float twodstatf::getStd2()
	{
		return quad/count-(sum1*sum1+sum2*sum2)/count/count;
	}
	float twodstatf::getMeanX()
	{
		return sum1/count;
	}
	float twodstatf::getMeanY()
	{
		return sum2/count;
	}
	long int twodstatf::getCount()
	{
		return count;
	}
	void twodstatf::add(float x, float y)
	{
		sum1 += x;
		sum2 += y;
		count++; 
		quad += x*x + y*y;
	}
	float twodstatf::getSum1()
	{
		return(sum1);
	}
	float twodstatf::getSum2()
	{
		return(sum2);
	}
	float twodstatf::getQuad()
	{
		return(quad);
	}


	//-------------------------------------------------------------------- stDev
	stDev::stDev()
	{
		reset();
	}
	stDev::~stDev()
	{
	}
	void stDev::set(stDev *devval)
	{
		quad = devval->getQuad(); 
		sum = devval->getSum(); 
		count = devval->getCount();
	}
	void stDev::add(stDev *devval)
	{
		quad += devval->getQuad(); 
		sum += devval->getSum(); 
		count += devval->getCount();
	}
	void stDev::reset()
	{
		count = 0; sum = 0.; quad = 0.;
	}
	double stDev::get()
	{
		if (count > 1)
		{
			double stdev2 = quad/(count-1)-sum*sum/count/(count-1);
			if (stdev2 > 0.)
				return sqrt(stdev2);
			else
				return 0.;
		}
		else
			return 0.;
	}
	double stDev::getErrorOfMean()
	{
		return get()/std::sqrt(double(count));
	}
	long int stDev::getCount()
	{
		return(count);
	}
	void stDev::add(double x, unsigned int n)
	{
		sum += x*n; 
		count += n; 
		quad += x*x*n;
	}
	double stDev::getQuad()
	{
		return(quad);
	}
	double stDev::getSum()
	{
		return(sum);
	}


	//--------------------------------------------------------------- C_Correlator
	C_Correlator::C_Correlator()
	{
		reset();
	}
	//copy constructor must be used, when using a vector<C_Correlator>
	//everything use with vector must have _type semantics_
	//i.e. return the exact same object when copied
	C_Correlator::C_Correlator( const C_Correlator &copy )
	{
		count = copy.count; 
		sum = copy.sum; 
		quad = copy.quad;
	}
	C_Correlator::~C_Correlator()
	{
	}
	
	//assignment operator, cf. above
	C_Correlator& C_Correlator::operator=(const C_Correlator &copy)
	{
		count = copy.count; 
		sum = copy.sum; 
		quad = copy.quad;
		
		return *this;
	}
	
	void C_Correlator::set(C_Correlator *devval)
	{
		quad = devval->getQuad(); 
		sum = devval->getSum(); 
		count = devval->getCount();
	}
	void C_Correlator::add(C_Correlator *devval)
	{
		quad += devval->getQuad(); 
		sum += devval->getSum(); 
		count += devval->getCount();
	}
	void C_Correlator::reset()
	{
		count = 0; 
		sum = 0.; 
		quad = 0.;
	}
	double C_Correlator::get()
	{
		return quad/count-sum*sum/count/count;
	}
	long int C_Correlator::getCount()
	{
		return count;
	}
	void C_Correlator::add(double x1, double x2)
	{
		sum += (x1+x2)/2.; 
		count++; 
		quad += x1*x2;
	}
	double C_Correlator::getQuad()
	{
		return quad;
	}
	double C_Correlator::getSum()
	{
		return sum;
	}


	//------------------------------------------------------------------- stDev2
	stdev2::stdev2()
	{
		reset();
	}
	stdev2::~stdev2()
	{
	}
	void stdev2::set(stdev2 *devval)
	{
		quad = devval->getQuad(); sum = devval->getSum(); 
		count = devval->getCount();
	}
	void stdev2::add(stdev2 *devval)
	{
		quad += devval->getQuad(); sum += devval->getSum(); 
		count += devval->getCount();
	}
	void stdev2::reset()
	{
		count = 0; sum = 0.; quad = 0.;
	}
	double stdev2::get()
	{
		return quad/count-sum*sum/count/count;
	}
	long int stdev2::getCount()
	{
		return count;
	}
	void stdev2::add(double x)
	{
		sum += x; count++; quad += x*x;
	}
	double stdev2::getQuad()
	{
		return quad;
	}
	double stdev2::getSum()
	{
		return sum;
	}


	//-------------------------------------------------------------------- statS
	statS::statS()
	{
		reset();
	}
	statS::~statS()
	{
	}
	int statS::reset()
	{
		count = 0;
		S = 0.;
		return 0;
	}
	long int statS::getCount()
	{
		return count;
	}
	int statS::add()
	{
		count++; S += 1./SIGMA/SIGMA;
		return 0;
	}
	int statS::add(double sig)
	{
		count++; S += 1./sig/sig;
		return(0);
	}
	double statS::get ()
	{
		return(S);
	}

	//------------------------------------------------------------------- statSx
	statSx::statSx()
	{
		reset();
	}
	statSx::~statSx()
	{
	}
	int statSx::reset()
	{
		count = 0; S = 0.;
		return 0;
	}
	long int statSx::getCount()
	{
		return count;
	}
	int statSx::add(double x)
	{
		count++; S += x/SIGMA/SIGMA;
		return 0;
	}
	int statSx::add(double x, double sig)
	{
		count++; S += x/sig/sig;
		return 0;
	}
	double statSx::get()
	{
		return(S);
	}


	//------------------------------------------------------------------- statSy
	statSy::statSy()
	{
		reset();
	}
	statSy::~statSy()
	{
	}
	int statSy::reset()
	{
		count = 0; S = 0.;
		return 0;
	}
	long int statSy::getCount()
	{
		return count;
	}
	int statSy::add(double y)
	{
		count++; S += y/SIGMA/SIGMA;
		return 0;
	}
	int statSy::add(double y, double sig)
	{
		count++;
		S += y/sig/sig;
		return 0;
	}
	double statSy::get ()
	{
		return S;
	}


	//------------------------------------------------------------------ statSxx
	statSxx::statSxx()
	{
		reset();
	}
	statSxx::~statSxx()
	{
	}
	int statSxx::reset()
	{
		count = 0;
		S = 0.;
		return 0;
	}
	long int statSxx::getCount()
	{
		return count;
	}
	int statSxx::add(double x)
	{
		count++;
		S += x*x/SIGMA/SIGMA;
		return 0;
	}
	int statSxx::add(double x, double sig)
	{
		count++; S += x*x/sig/sig;
		return 0;
	}
	double statSxx::get()
	{
		return S;
	}

	//------------------------------------------------------------------ statSxy
	statSxy::statSxy()
	{
		reset();
	}
	statSxy::~statSxy()
	{
	}
	int statSxy::reset()
	{
		count = 0;
		S = 0.;
		return 0;
	}
	long int statSxy::getCount()
	{
		return count;
	}
	int statSxy::add(double x, double y)
	{
		count++;
		S +=x*y/SIGMA/SIGMA;
		return 0;
	}
	int statSxy::add(double x, double y, double sig)
	{
		count++;
		S +=x*y/sig/sig;
		return 0;
	}
	double statSxy::get ()
	{
		return(S);
	}


	//------------------------------------------------------------------ linFit
	linFit::linFit()
	{
		S.reset();
		Sx.reset();
		Sy.reset();
		Sxx.reset();
		Sxy.reset();
		Delta = 0.;
		a = 0.;
		b = 0.;
	}
	int linFit::add(double xval, double yval)
	{
		S.add();
		Sx.add(xval);
		Sy.add(yval);
		Sxx.add(xval);	
		Sxy.add(xval, yval);
		return 0;
	}
	int linFit::add(double xval, double yval, double sig)
	{
		S.add(sig);
		Sx.add(xval, sig);
		Sy.add(yval, sig);
		Sxx.add(xval, sig);	
		Sxy.add(xval, yval, sig);
		return 0;
	}
	double linFit::getA()
	{
		Delta = S.get()*Sxx.get() - Sx.get()*Sx.get();
		a = (Sxx.get()*Sy.get() - Sx.get()*Sxy.get())/Delta;
		return a;
	}
	double linFit::getB()
	{
		Delta = S.get()*Sxx.get() - Sx.get()*Sx.get();
		b = (S.get()*Sxy.get() - Sx.get()*Sy.get())/Delta;
		return b;
	}
	double linFit::getSlope()
	{
		return getB();
	}
	double linFit::getOffset()
	{
		return getA();
	}


	//------------------------------------------------------------------ Queue
	Queue::Queue(long int lengthval)
	{
    	length = lengthval;
    	counter = 0;
    	value = new double[length];
    	for (int i = 0; i < length; i++)
		value[i] = 0.;
	}

	Queue::Queue(long int lengthval, double initval)
	{
    	length = lengthval;
    	counter = 0;
    	value = new double[length];
    	for (int i = 0; i < length; i++)
		value[i] = initval;
	}

	Queue::~Queue()
	{
    	delete []value;
	}

	long int Queue::getLength()
	{
    	return length;
	}

	long int Queue::getCount()
	{
    	return counter;
	}

	double Queue::getVal()
	{
    	return getVal(length-1);
	}

	double Queue::getVal(long int i)
	{
    	return value[i];
	}

	void Queue::add(double val)
	{
    	for (long int i = length-2; i >= 0; i--)
		{
			value[i+1] = value[i];
		}
    	counter++;
    	value[0] = val;
	}

	double Queue::getMean()
	{
    	double sum = 0.; 
    	for (int i = 0; (i < length) && (i < counter); i++)
		{
			sum += value[i];
		}
    	sum /= (length > counter ? counter : length);
    	return sum;
	}

	double Queue::getDiff()
	{
    	double diff = value[(length > counter ? counter-1 : length - 1)] - value[0];
    	return diff;
	}

	double Queue::getLast()
	{
    	return getVal(0);
	}


	//==========================================================================
	//                                                             interpolation
	//==========================================================================
	
	double interpolate(double val0, double val1, double t)
	//	interpolates linearly between the values at 0 and at 1 using the 
	//	parameter 0 <= t <= 1 to indicate the position of the interpolated 
	//	point between 0 and 1. 
	//	If t lies out of bounds, it is silently set to the closest position
	{
		if (t < 0.) t = 0.;
		if (t > 1.) t = 1.;
		return ((1.-t)*val0 + t*val1);
	}
	
	double interpolate2D(double val00, double val01, double val10, double val11, 
						 double ty, double tx)
	//	interpolates linearly in 2D. The first index runs in y, the second in x. 
	//	ty and tx indicate the postion of the interpolation point in the two 
	//	dimensions
	{
		double val0 = interpolate(val00, val01, tx);
		double val1 = interpolate(val10, val11, tx);
		return(interpolate(val0, val1, ty));
	}
	
	void absphase2reim(double abs, double phase, double &re, double &im)
	{
		re = abs*cos(phase);
		im = abs*sin(phase);
	}
	
	void reim2absphase(double re, double im, double &abs, double &phase)
	{
		abs = sqrt(re*re + im*im);
		if (re == 0.)
		{
			if (im == 0.)
				phase = 0.;
			else 
			{
				if (im > 0.)
					phase = pi/2.;
				else
					phase = -pi/2.;
			}
		}
		
		if (re > 0.) // quadrant 2 and 4
			phase = atan(im/re);
		if (re < 0.) 
		{
			if (im > 0.) // quadrant 1
				phase = pi + atan(im/re);
			else if(im < 0.) // quadrant 3
				phase = -pi + atan(im/re);
		}
	}

	void interpolate(double valre0, double valim0, double valre1, double valim1, double t, double &resre, double &resim)
	{
		resre = interpolate(valre0, valre1, t);
		resim = interpolate(valim0, valim1, t);
	}
	
	void interpolatePhase(double valre0, double valim0, double valre1, double valim1, double t, double &resre, double &resim)
	{
		double tmpre = interpolate(valre0, valre1, t);
		double tmpim = interpolate(valim0, valim1, t);

		double abs0 = sqrt(valre0*valre0 + valim0*valim0);
		double abs1 = sqrt(valre1*valre1 + valim1*valim1);

		double abs = interpolate(abs0, abs1, t);		
		double abstmp = sqrt(tmpre*tmpre + tmpim*tmpim);
		
		resre = abs/abstmp*tmpre;
		resim = abs/abstmp*tmpim;
	}
	
	void interpolate2D(double valre00, double valim00, double valre01, double valim01, double valre10, double valim10, 
					   double valre11, double valim11, double ty, double tx, double &resre, double &resim)
	{
		double resre0, resim0, resre1, resim1;
		interpolate(valre00, valim00, valre01, valim01, tx, resre0, resim0);
		interpolate(valre10, valim10, valre11, valim11, tx, resre1, resim1);
		interpolate(resre0, resim0, resre1, resim1, ty, resre, resim);
	}
	
	void interpolate2DPhase(double valre00, double valim00, double valre01, double valim01, double valre10, double valim10, 
					   double valre11, double valim11, double ty, double tx, double &resre, double &resim)
	{
		double resre0, resim0, resre1, resim1;
		interpolatePhase(valre00, valim00, valre01, valim01, tx, resre0, resim0);
		interpolatePhase(valre10, valim10, valre11, valim11, tx, resre1, resim1);
		interpolatePhase(resre0, resim0, resre1, resim1, ty, resre, resim);
	}
	
	//==========================================================================
	//                                                             miscellaneous
	//==========================================================================

//***commented out the following class for the use in the cross correlator project, but matlib generally needs this
/*
	//--------------------------------------------------------------- C_Progress
	C_Progress::C_Progress()
	{
		init();
	}

	C_Progress::C_Progress( double max_value )
	{
		init( max_value );
	}
	
	C_Progress::C_Progress( std::istream &in )
	{
		init( in );
	}

	C_Progress::C_Progress( std::ostream &out )
	{
		init( out );
	}

	C_Progress::C_Progress( const C_Progress &that )
	{
		p_progress_0	= that.p_progress_0;
		p_progress_1	= that.p_progress_1;
		p_progress		= that.p_progress;
		p_timestep		= that.p_timestep;
		p_prefix		= that.p_prefix;
		p_timer 		= that.p_timer;
	}
	
	C_Progress::~C_Progress()
	{
	}

	C_Progress &C_Progress::operator=( const C_Progress &that )
	{
		p_progress_0	= that.p_progress_0;
		p_progress_1	= that.p_progress_1;
		p_progress		= that.p_progress;
		p_timestep		= that.p_timestep;
		p_prefix		= that.p_prefix;
		p_timer 		= that.p_timer;
		return *this;
	}
	
	void C_Progress::set_timestep( double timestep )
	{
		p_timestep = timestep;
	}
	
	void C_Progress::set_prefix( string prefix )
	{
		p_prefix = prefix;
	}
	
	double C_Progress::get_timestep()
	{
		return p_timestep;
	}
	
	string C_Progress::get_prefix()
	{
		return p_prefix;
	}
	
	void C_Progress::init()
	{
		p_progress_0 = 0;	// Not needed; just for aesthetic reasons...
		p_progress_1 = 0;	// Not needed; just for aesthetic reasons...
		p_progress = 0;
		p_prefix = DEFAULT_PREFIX;
		p_timestep = DEFAULT_TIMESTEP;
		p_timer.set_timer( get_timestep(), ns_timer::msec );
	}
	
	void C_Progress::init( double max_value )
	{
		p_progress_0 = 0;
		p_progress_1 = max_value;
		p_prefix = DEFAULT_PREFIX;
		p_timestep = DEFAULT_TIMESTEP;
		p_timer.set_timer( get_timestep(), ns_timer::msec );
		ns_commandline::cursorOff();
		cout << get_prefix() << "   0%" << std::flush;
	}
	
	void C_Progress::init( std::istream &in )
	{
		std::ios::pos_type pos0;
		
		pos0 = in.tellg();
		p_progress_0 = pos0;
		in.seekg( 0, std::ios::end );
		p_progress_1 = in.tellg();
		in.seekg( pos0 );
		ns_commandline::cursorOff();
		p_prefix = DEFAULT_PREFIX;
		p_timestep = DEFAULT_TIMESTEP;
		p_timer.set_timer( get_timestep(), ns_timer::msec );
		cout << get_prefix() << "   0%" << std::flush;
	}

	void C_Progress::init( std::ostream &out )
	{
		std::ios::pos_type pos0;
		
		pos0 = out.tellp();
		p_progress_0 = pos0;
		out.seekp( 0, std::ios::end );
		p_progress_1 = out.tellp();
		out.seekp( pos0 );
		ns_commandline::cursorOff();
		p_prefix = DEFAULT_PREFIX;
		p_timestep = DEFAULT_TIMESTEP;
		p_timer.set_timer( get_timestep(), ns_timer::msec );
		cout << get_prefix() << "   0%" << std::flush;
	}

	void C_Progress::show( int units_per_asterisk, int asterisks_per_line )
	{
		int units_per_line;
		
		units_per_line = units_per_asterisk*asterisks_per_line;
	
		if ( p_progress % units_per_line == 0 )
		{
			cout << endl;
			cout << "Time: " << clock()/CLOCKS_PER_SEC << "s\t";
			cout << p_progress << "\t" << std::flush;
		}

		p_progress++;

		if ( p_progress % units_per_asterisk == 0 )
		{
			cout << "*" << std::flush;
		}
	}

	void C_Progress::show( double value )
	{
		double progress;
		double percent;

		if ( p_timer.timeout() )
		{
			progress = value;
			percent = ceil( 100*double(progress-p_progress_0)/double(p_progress_1-p_progress_0) );

			cout << char(8) << char(8) << char(8) << char(8) << char(8) << std::fixed << std::setprecision(0) << std::setw( 4 )
				 << percent << std::setw(0) << "%" << std::flush;
			cout.unsetf(std::ios_base::fixed | std::ios_base::scientific);
			cout << std::setprecision(6);
			p_timer.set_timer( get_timestep(), ns_timer::msec );
		}
	}
	
	void C_Progress::show( std::istream &in )
	{
		double progress;
		double percent;

		if ( p_timer.timeout() )
		{
			progress = in.tellg();
			percent = ceil( 100*double(progress-p_progress_0)/double(p_progress_1-p_progress_0) );

			cout << char(8) << char(8) << char(8) << char(8) << char(8) << std::fixed << std::setprecision(0) << std::setw( 4 )
				 << percent << std::setw(0) << "%" << std::flush;
			cout.unsetf(std::ios_base::fixed | std::ios_base::scientific);
			cout << std::setprecision(6);

			p_timer.set_timer( get_timestep(), ns_timer::msec );

		}
	}
	
	void C_Progress::show( std::ostream &out )
	{
		double progress;
		double percent;
		
		if ( p_timer.timeout() )
		{
			progress = out.tellp();
			percent = ceil( 100*double(progress-p_progress_0)/double(p_progress_1-p_progress_0) );

			cout << char(8) << char(8) << char(8) << char(8) << char(8) << std::fixed << std::setprecision(0) << std::setw( 4 )
				 << percent << std::setw(0) << "%" << std::flush;
			cout.unsetf(std::ios_base::fixed | std::ios_base::scientific);
			cout << std::setprecision(6);

			p_timer.set_timer( get_timestep(), ns_timer::msec );

		}
	}
	
	void C_Progress::finished()
	{
		cout << char(8) << char(8) << char(8) << char(8) << "100%" << endl;
		ns_commandline::cursorOn();
	}
*/


	//--------------------------------------------------------- C_LoopDataSingle
	C_LoopDataSingle::C_LoopDataSingle()
	{
		init();
	}
	
	C_LoopDataSingle::C_LoopDataSingle( double start, double stop, int steps )
	{
		init();
		set( start, stop, steps );
	}
	
	C_LoopDataSingle::C_LoopDataSingle( double start, double stop, int steps, loopfunc_t func )
	{
		init();
		set( start, stop, steps, func );
	}
	
	C_LoopDataSingle::C_LoopDataSingle( const C_LoopDataSingle &loopdatasingle )
	{
		p_start = loopdatasingle.p_start;
		p_stop = loopdatasingle.p_stop;
		p_steps = loopdatasingle.p_steps;
		p_count = loopdatasingle.p_count;
		p_func = loopdatasingle.p_func;
		p_values = loopdatasingle.p_values;
	}
	
	C_LoopDataSingle::~C_LoopDataSingle()
	{
	}
	
	C_LoopDataSingle &C_LoopDataSingle::operator=( const C_LoopDataSingle &loopdatasingle )
	{
		p_start = loopdatasingle.p_start;
		p_stop = loopdatasingle.p_stop;
		p_steps = loopdatasingle.p_steps;
		p_count = loopdatasingle.p_count;
		p_func = loopdatasingle.p_func;
		p_values = loopdatasingle.p_values;
		return *this;
	}
	
	void C_LoopDataSingle::init()
	{
		p_start = 0;
		p_stop = 0;
		p_steps = 0;
		p_count = 0;
		p_func = 0;
		p_values.clear();
		p_values.push_back( get_value() );
	}

	void C_LoopDataSingle::set_start( double start )
	{
		p_start = start;
		reset();
	}
	
	void C_LoopDataSingle::set_stop( double stop )
	{
		p_stop = stop;
		reset();
	}
	
	void C_LoopDataSingle::set_steps( int steps )
	{
		if ( steps < 0 )
		{
			p_steps = 0;
		}
		else
		{
			p_steps = steps;
		}
		reset();
	}

	void C_LoopDataSingle::set_func( loopfunc_t func )
	{
		p_func = func;
		reset();
	}
	
	void C_LoopDataSingle::set( double start, double stop, int steps )
	{
		set_start( start );
		set_stop( stop );
		set_steps( steps );
	}

	void C_LoopDataSingle::set( double start, double stop, int steps, loopfunc_t func )
	{
		set_start( start );
		set_stop( stop );
		set_steps( steps );
		set_func( func );
	}
	
	void C_LoopDataSingle::reset()
	{
		p_count = 0;
		p_values.clear();
		p_values.push_back( get_value() );
	}

	double C_LoopDataSingle::get_start() const
	{
		return p_start;
	}
	
	double C_LoopDataSingle::get_stop() const
	{
		return p_stop;
	}
	
	int C_LoopDataSingle::get_steps() const
	{
		return p_steps;
	}

	loopfunc_t C_LoopDataSingle::get_func() const
	{
		return p_func;
	}

	bool C_LoopDataSingle::next()
	{
		p_count++;
		if ( p_count <= p_steps )
		{
			p_values.push_back( get_value() );
			return true;
		}
		else
		{
			return false;
		}
	}
	
	double C_LoopDataSingle::get_value() const
	{
		double value;
		
		if ( get_steps() != 0 )
		{
			value = get_start() + ( get_stop() - get_start() )/get_steps() * get_count();
		}
		else
		{
			value = get_start();
		}
		
		if ( p_func != 0 )
		{
			value = p_func( value );
		}
		
		return value;
	}
	
	double C_LoopDataSingle::get_sum( int depth ) const
	{
		double sum;

		sum = 0;
		for ( int i=0; i<(int )p_values.size()-depth; i++ )
		{
			sum += p_values[i];
		}		

		return sum;
	}
	
	double C_LoopDataSingle::get_sum_range( int start, int stop )
	{
		double sum=0;

		cout << "p_steps=" << p_steps << endl;
		
		reset();
		do
		{
			cout << "--> p_values.size()=" << (int)p_values.size() << endl;
		} while ( next() );
		cout << "==> p_values.size()=" << (int)p_values.size() << endl;
		
		if ( start < 0 )
		{
			start = 0;
		}
		
		if ( ( stop >= (int )p_values.size() ) || ( stop < 0 ) )
		{
			stop = (int)p_values.size() - 1;
		}
		
		sum = 0;
		for ( int i=start; i<=stop; i++ )
		{
			sum += p_values[i];
		}

		cout << "p_values.size()=" << (int)p_values.size() << endl;

		return sum;
	}

	int C_LoopDataSingle::get_count() const
	{
		return p_count;
	}
	


	//--------------------------------------------------------------- C_LoopData
	C_LoopData::C_LoopData()
	{
		init();
	}
	
	C_LoopData::C_LoopData( const C_LoopData &loopdata )
	{
		p_looplist = loopdata.p_looplist;
		p_count = loopdata.p_count;
	}
	
	C_LoopData::~C_LoopData()
	{
	}

	C_LoopData &C_LoopData::operator=( const C_LoopData &loopdata )
	{
		p_looplist = loopdata.p_looplist;
		p_count = loopdata.p_count;
		return *this;
	}

	void C_LoopData::init()
	{
		p_looplist.clear();
		p_count = 0;
	}

	void C_LoopData::clear()
	{
		p_looplist.clear();
		reset();
	}

	void C_LoopData::reverse()
	{
		std::reverse( p_looplist.begin(), p_looplist.end() );
		reset();
	}
	
	void C_LoopData::add( const C_LoopDataSingle &loopdatasingle )
	{
		p_looplist.push_back( loopdatasingle );
		reset();
	}
	
	void C_LoopData::add( double start, double stop, int steps )
	{
		add( C_LoopDataSingle( start, stop, steps ) );
	}
	
	void C_LoopData::add( double start, double stop, int steps, loopfunc_t func )
	{
		add( C_LoopDataSingle( start, stop, steps, func ) );
	}
	
	void C_LoopData::add_loops( looplist_t &looplist )
	{
		p_looplist.insert( p_looplist.end(), looplist.begin(), looplist.end() );
		reset();
	}
	
	C_LoopDataSingle C_LoopData::get_loop_copy( int index ) const
	{
		C_LoopDataSingle loop;
		
		if ( index < get_dimension() )
		{
			loop = p_looplist[index];
		}
		else
		{
			loop.init();
		}
		
		return loop;
	}

	C_LoopDataSingle &C_LoopData::get_loop( int index )
	{
		if ( index < get_dimension() )
		{
			return p_looplist[index];
		}
		else
		{
			C_LoopDataSingle *pntr = new C_LoopDataSingle();
			return *pntr;
		}
	}

	double C_LoopData::get_start( int index ) const
	{
		double start;
		
		if ( index < get_dimension() )
		{
			start = p_looplist[index].get_start();
		}
		else
		{
			start = 0;
		}
		
		return start;
	}
	
	double C_LoopData::get_stop( int index ) const
	{
		double stop;
		
		if ( index < get_dimension() )
		{
			stop = p_looplist[index].get_stop();
		}
		else
		{
			stop = 0;
		}
		
		return stop;
	}
	
	int C_LoopData::get_steps( int index ) const
	{
		int steps;
		
		if ( index < get_dimension() )
		{
			steps = p_looplist[index].get_steps();
		}
		else
		{
			steps = 0;
		}
		
		return steps;
	}
	
	loopfunc_t C_LoopData::get_func( int index ) const
	{
		loopfunc_t func;
		
		if ( index < get_dimension() )
		{
			func = p_looplist[index].get_func();
		}
		else
		{
			func = 0;
		}
		
		return func;
	}
	
	double C_LoopData::get_value( int index ) const
	{
		double value;
		
		if ( index < get_dimension() )
		{
			value = p_looplist[index].get_value();
		}
		else
		{
			value = 0;
		}
		
		return value;
	}
	
	double C_LoopData::get_sum( int index, int depth ) const
	{
		double sum;
		
		if ( index < get_dimension() )
		{
			sum = p_looplist[index].get_sum( depth );
		}
		else
		{
			sum = 0;
		}
		
		return sum;
	}
	
	double C_LoopData::get_sum_range( int index, int start, int stop )
	{
		double sum;
		
		if ( index < get_dimension() )
		{
			sum = p_looplist[index].get_sum_range( start, stop );
		}
		else
		{
			sum = 0;
		}
		
		return sum;
	}
	
	int C_LoopData::get_count( int index ) const
	{
		int count;
		
		if ( index < get_dimension() )
		{
			count = p_looplist[index].get_count();
		}
		else
		{
			count = 0;
		}
		
		return count;

	}

	int C_LoopData::get_count() const
	{
		return p_count;
	}

	int C_LoopData::get_dimension() const
	{
		return (int)p_looplist.size();
	}

	bool C_LoopData::next()
	{
		p_count++;
		
		for ( int i=get_dimension()-1; i>=0; i-- )
		{
			if ( p_looplist[i].next() )
			{
				return true;
			}
			p_looplist[i].reset();
		}
		
		return false;
	}
	
	void C_LoopData::reset()
	{
		for ( int i=0; i<get_dimension(); i++ )
		{
			p_looplist[i].reset();
		}
		p_count = 0;
	}

} // namespace
