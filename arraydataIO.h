//
//  arraydataIO.h
//  xcca_commandline
//
//  Created by Feldkamp on 7/3/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//
//
//	handles input/output to/from arraydata objects
//
//	this class relies on other libraries to provide the io for different formats
//	the dependence can be controlled by the ARRAYDATAIO_xxx preprocessor macros 
//	if none of the i/o is actually needed, you can undefine all of the macros
	#define ARRAYDATAIO_EDF				//use EDF
	#define ARRAYDATAIO_TIFF			//use TIFF
	#define ARRAYDATAIO_HDF5			//use HDF5

//	#undef ARRAYDATAIO_EDF				//do not use EDF
//	#undef ARRAYDATAIO_TIFF				//do not use TIFF
//	#undef ARRAYDATAIO_HDF5				//do not use HDF5

//  array2D<double> assumes matrix convention, i.e. dimensions (rows, columns).
//  The read/write operations assume image convention (nx, ny).
//  The transpose flag in each function can be used to
//  transpose the image before reading or after writing
//  --- rule of thumb: generally, transpose should be true ---

#ifndef _arraydataIO_h
#define _arraydataIO_h

	#include <string>
	#include "arrayclasses.h"
	using namespace ns_arraydata;

	class arraydataIO{

	public:
		arraydataIO( int verbose = 0 );
		~arraydataIO();
		
		//general functions, data type is determined by checking the file extension
		int readFromFile( std::string filename, array1D<double> *&dest) const;
		int readFromFile( std::string filename, array2D<double> *&dest) const;
		int writeToFile( std::string filename, array1D<double> *src) const;
		int writeToFile( std::string filename, array2D<double> *src) const;

		//------------------------------------------------------------------------------EDF
		int readFromEDF( std::string filename, array1D<double> *&dest ) const;
		int readFromEDF( std::string filename, array2D<double> *&dest ) const;
		int writeToEDF( std::string filename, array1D<double> *src ) const;
		int writeToEDF( std::string filename, array2D<double> *src ) const;
		
		//------------------------------------------------------------------------------Tiff	
		int readFromTiff( std::string filename, array2D<double> *&dest ) const;		// there is a weird problem with reading unscaled data, scaled is fine
		
		// scaleFlag = 0 --> do not scale data
		// scaleFlag = 1 --> scale data to 0..65535
		int writeToTiff( std::string filename, array2D<double> *src, int scaleFlag = 0, int bits = 16) const;
		
		
		//------------------------------------------------------------------------------HDF5
		int readFromHDF5( std::string filename, array1D<double> *&dest) const;
		int readFromHDF5( std::string filename, array2D<double> *&dest) const;
		
		// dataType = 0 --> write as doubles (H5T_NATIVE_DOUBLE)
		// dataType = 1 --> write as float   (H5T_NATIVE_FLOAT, default)
		// dataType = 2 --> write as int     (H5T_NATIVE_INT)
		// dataType = 3 --> write as int16_t (H5T_STD_I16LE)
		// dataType = 4 --> write as long    (H5T_NATIVE_LONG)
		int writeToHDF5( std::string filename, array1D<double> *src, int dataType = 1) const;
		int writeToHDF5( std::string filename, array2D<double> *src, int dataType = 1) const;	
		
		//------------------------------------------------------------------------------ASCII
		int readFromASCII( std::string filename, array1D<double> *&dest ) const;
		int readFromASCII( std::string filename, array2D<double> *&dest ) const;
		int writeToASCII( std::string filename, array1D<double> *src, int format=0 ) const;
		int writeToASCII( std::string filename, array2D<double> *src, int format=0 ) const;
		int writeToASCII( std::string filename, array3D<double> *src, int format=0 ) const;


		//------------------------------------------------------------------------------general support functions
		bool transposeForIO() const;
		void setTransposeForIO( bool t );
		int verbose() const;
		void setVerbose( int level );
		
		void generateTestPattern( array2D<double> *img, int type );
		
	private:
		bool p_transposeForIO;	 // transpose 2D files to conform with matrix (rows, cols) convention (default: true)
		int p_verbose;
		
	};
	
	


#endif
