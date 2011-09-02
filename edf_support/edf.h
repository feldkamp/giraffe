//==========================================================================
//	Title: edf.h
//
//	Author: Christian G. Schroer
//
//	description: Class for ESRF-Data Format 
//
//--------------------------------------------------------------------------
// 
// http://www.esrf.eu/computing/expg/subgroups/general/format/Format.html
//
// For a list of predefined header keywords see file <edfkeywords.h>
//
//==========================================================================

#ifndef edf_h_flag
#define edf_h_flag

#include <string>
#include <sstream>
#include <vector>

#include "edfkeywords.h"
#include "timer.h"

#ifndef tomo_version
#define tomo_version "tomo_version_9.0"
#endif

namespace ns_edf
{
	//------------------------------------------------------------------------------
	// The following string constant should be always increased when manipulating
	// one of the edf source files edf.cpp, edf.h, or edfkeywords.h.
	//
	//const char EDF_VERSION_NUMBER[] = "edf_version_10.0";
	//
	// (by JP)
	// Bug removal: difference between 32bit and 64bit machines.
	//const char EDF_VERSION_NUMBER[] = "edf_version_10.1"; 
	//
	// (2009-09-02, JP) Implementation of complex data types and (c)volumedata.
	//const char EDF_VERSION_NUMBER[] = "edf_version_10.2"; 
	//
	// (2011-05-19, JP)
	// Bug removal: Introduction of matlibFlag in order to distinguish between
	// old matlib edf files and non-matlib edf files which is relevant for
	// backward compatibility (errornous byte order handling of old matlib edf
	// files).
	const char EDF_VERSION_NUMBER[] = "edf_version_10.3";
	
	
	//------------------------------------------------------------------------------
	// Constants used for parsing and writing edf header
	//
	const int HEADER_WIDTH_KEYWORD			= 20;					// Minimal number of characters for keyword.
	const int HEADER_WIDTH_VALUE 			= 23;					// Minimal number of characters for value.
	const char HEADER_TOKEN_BEGIN[] 		= "{";					// Token indicating begin of header.
	const char HEADER_TOKEN_END[] 			= "}";					// Token indicating end of header.
	const char HEADER_TOKEN_ASSIGN[]  	 	= " = ";    			// Token between tag and value.
	const char HEADER_TOKEN_SEPARATOR[]	  	= " ;";	  				// Token indicating begin of comments
	const char HEADER_TOKEN_SEPARATOR_ALT[]	= ";";	  				// Alternative token indicating begin of comments
																	// (only used for parsing, but not for writing).
	const std::string HEADER_TOKEN_ENDLINE	= std::string( 1, 10 ); // Character indicating end of header line.
	const int HEADER_BYTES_MODULO			= 1024; 				// Header size must be a multiple of this number.
	
	//------------------------------------------------------------------------------
	// Constants used for parsing and writing header ID.
	//
	const char HEADER_ID_PREFIX[]			= "EH";
	const char HEADER_ID_SEPARATOR[]		= ":";
	const int HEADER_ID_FIELDSIZE			= 6;
	
	
	//------------------------------------------------------------------------------
	// Special string values for header.
	//
	const char NOT_SPECIFIED_STR_TUD[]	= "NotSpecified";		// Property not defined (old)
	const char NOT_SPECIFIED_STR_ESRF[]	= "NoSpecificValue";	// Property not defined (new as given by edf documentation)
	const char FAMELESS_STR[]			= "Fameless";			// Property defined but invalid value.
	
	        
	//------------------------------------------------------------------------------
	// Byte order string values for header.
	//
	const char HIGH_BYTE_FIRST_STR[]	= "HighByteFirst";		// Most significant byte at lowest address (big endian).
	const char LOW_BYTE_FIRST_STR[]		= "LowByteFirst";		// Less significant byte at lowest address (little endian).
	
	
	//------------------------------------------------------------------------------
	// Data type string values for header.
	//
	const char FLOAT_SINGLE_STR[]		= "FloatValue"; 			// Floating point single precision (4 Bytes).
	const char FLOAT_DOUBLE_STR[] 		= "DoubleValue";			// Floating point single precision (8 Bytes).
	const char BYTE_STR[]				= "SignedByte"; 			// Signed char (1 Byte).
	const char SHORT_STR[]				= "SignedShort";			// Signed short (2 Bytes).
	const char INTEGER_STR[]			= "SignedInteger";			// Signed integer (4 Bytes).
	const char LONG_STR[]				= "SignedLong"; 			// Signed long (4 Bytes).
	const char UNSIGNED_BYTE_STR[]		= "UnsignedByte";			// Unsigned char (1 Byte).
	const char UNSIGNED_SHORT_STR[] 	= "UnsignedShort";			// Unsigned short (2 Bytes).
	const char UNSIGNED_INTEGER_STR[]	= "UnsignedInteger";		// Unsigned integer (4 Bytes).
	const char UNSIGNED_LONG_STR[]		= "UnsignedLong";			// Unsigned long (4 Bytes).
	const char COMPLEX_SINGLE_STR[]		= "ComplexFloat";			// 2 x float (8 Bytes)
	const char COMPLEX_DOUBLE_STR[] 	= "ComplexDouble";			// 2 x double (16 Bytes)
	const char COMPLEX_BYTE_STR[]		= "ComplexByte";			// 2 X signed byte (2 Bytes)
	const char COMPLEX_SHORT_STR[]		= "ComplexShortInteger";	// 2 x signed short (4 Bytes)
	const char COMPLEX_INTEGER_STR[]	= "ComplexInteger"; 		// 2 x signed integer (8 Bytes)
	const char COMPLEX_LONG_STR[]		= "ComplexLongInteger"; 	// 2 x signed long (8 Bytes)
	
	//------------------------------------------------------------------------------
	// Deprecated data type string values kept for backwards compatibility.
	//
	const char OLD_FLOAT_STR[]			= "Float";			// Same as FloatValue.
	const char OLD_DOUBLE_STR[]			= "Double"; 		// Same as DoubleValue.
	const char OLD_SHORT_STR[]			= "Short";			// Same as SignedShort.
	const char OLD_LONG_STR[] 			= "Long";			// Same as SignedInteger.
	const char OLD_UNSIGNED_SHORT_STR[] = "UnsignedShort";	// Same as UnsignedShort.
	const char OLD_UNSIGNED_LONG_STR[]	= "UnsignedLong";	// Same as UnsignedInteger.
	
	
	//------------------------------------------------------------------------------
	// Real format string values for header.
	//
	const char REAL_FORMAT_IEEE_STR[]		= "IEEE";		// IEEE 754 standard.
	const char REAL_FORMAT_VAX_STR[]		= "VAX";		// Not implemented.
	const char REAL_FORMAT_VMS_STR[]		= "VMS";		// Not implemented.
	const char REAL_FORMAT_CONVEX_STR[]	= "ConvexReal"; 	// Not implemented.
	

	//--------------------------------------------------------------------------
	// File type string values for header.
	//
	const char FTS_VECTOR_DATA[]		  = "vectordata";			//
	const char FTS_MATRIX_DATA[]		  = "matrixdata";			//
	const char FTS_VOLUME_DATA[]		  = "volumedata";			//
	const char FTS_CVECTOR_DATA[]		  = "cvectordata";			//
	const char FTS_CMATRIX_DATA[]		  = "cmatrixdata";			//
	const char FTS_CVOLUME_DATA[]		  = "cvolumedata";			//
	const char FTS_IMAGE_DATA[] 		  = "imagedata";			//
	const char FTS_SINO_DATA[]  		  = "sinodata"; 			//
	const char FTS_MU_DATA[]			  = "mudata";				//
	const char FTS_MATRIX_STACK[]		  = "matrix_stack";  		//
	const char FTS_IMAGE_STACK[]		  = "image_stack"; 		  	//
	const char FTS_SINO_STACK[] 		  = "sino_stack";  		  	//
	const char FTS_FLUO_TOMO_IMAGE[]	  = "fluotomo_image";		//
	const char FTS_ATTENUATION_SINOGRAM[] = "attenuation_sinogram"; //
	const char FTS_DIFFERENCE_SINOGRAM[]  = "difference_sinogram";  //
	const char FTS_GENERIC_EDF[]		  = "generic_edf";			//
	const char FTS_GENERIC_VOXEDF[] 	  = "generic_voxedf";  		//
	
	const char FTS_BINNED[]				  = "binned";
	const char FTS_BINNNED_EDF[]		  = "binned_edf";
	const char FTS_BINNED_SINOGRAM[]	  = "binned_sinogram";
	const char FTS_FBT_REC[]			  = "fbt_rec";
	const char FTS_FFBT_REC[] 			  = "ffbt_rec";
	const char FTS_FILTERED_SINOGRAM[]	  = "filtered_sinogram";
	const char FTS_FLIPPED_SINOGRAM[] 	  = "flipped_sinogram";
	const char FTS_FLIPPED[]			  = "flipped";
	const char FTS_HORIZ_NORMALIZED[] 	  = "horiz_normalized";
	const char FTS_IMAGE_SUPPORT[]		  = "image_support";
	const char FTS_POINTWISE_ABS[]		  = "pointwise_abs";
	const char FTS_RESCALE_REGION[]		  = "rescale_region";
	const char FTS_SHIFTED_SINO[] 		  = "shifted_sino";
	const char FTS_SINOGRAM[] 			  = "sinogram";
	const char FTS_SINOGRAM_ABS[] 		  = "sinogram_abs";
	const char FTS_SINOGRAM_FAN[] 		  = "sinogram_fan";
	const char FTS_SUB_MEAN_NOISE[] 	  = "sub_mean_noise";
	const char FTS_SINOGRAM_SUPPORT[] 	  = "sinogram_support";
	const char FTS_SINOGRAM_FF_BY_CONV[]  = "sinogram_ff_by_conv";
	const char FTS_THRESHOLDED[]  		  = "thresholded";
	const char FTS_ZOOMED[]				  = "zoomed";


	// ByteOrderOutMode
	enum boom_t
	{
		BOOM_INPUT, 	// Output byte order determined by input.
		BOOM_MACHINE,	// Output byte order determined by machine.
		BOOM_USER,		// Output byte order determined by user.
	};


	enum byteorder_t
	{
		BO_HIGH_BYTE_FIRST  =  0,
		BO_LOW_BYTE_FIRST	=  1,
		
		BO_NOT_SPECIFIED	= -1,
		BO_FAMELESS 		= -2,
	};


	enum datatype_t
	{
		DT_FLOAT			=  0,
		DT_DOUBLE			=  1,
		
		DT_BYTE				=  2,
		DT_SHORT			=  3,
		DT_INTEGER  		=  4,
		DT_LONG				=  5,
		
		DT_UNSIGNED_BYTE	=  6,
		DT_UNSIGNED_SHORT	=  7,
		DT_UNSIGNED_INTEGER	=  8,
		DT_UNSIGNED_LONG	=  9,
		
		DT_COMPLEX_FLOAT	= 10,
		DT_COMPLEX_DOUBLE	= 11,
		DT_COMPLEX_BYTE 	= 12,
		DT_COMPLEX_SHORT	= 13,
		DT_COMPLEX_INTEGER	= 14,
		DT_COMPLEX_LONG 	= 15,
		
		DT_NOT_SPECIFIED	= -1,
		DT_FAMELESS 		= -2,
	};
	const int datatype_min =  0; 	// Useful for iterating over datatype objects in a loop.
	const int datatype_max = 15; 	// Useful for iterating over datatype objects in a loop.
	std::string get_datatypes();
	std::string get_datatype_string( ns_edf::datatype_t type );


	enum realformat_t
	{
		REAL_FORMAT_IEEE	= 0,
		REAL_FORMAT_VAX 	= 1,
		REAL_FORMAT_VMS 	= 2,
		REAL_FORMAT_CONVEX	= 3,
		
		REAL_FORMAT_NOT_SPECIFIED = -1,
		REAL_FORMAT_FAMELESS	  = -2,
	};


	enum scaled_t
	{
		SF_UNSCALED = 0,
		SF_SCALED	= 1,
		SF_RETAIN	= 2,
	};


	typedef std::string filetype_t;
/*
	enum filetype_t
	{
		FT_VECTOR_DATA			=  0,
		FT_MATRIX_DATA			=  1,
		FT_IMAGE_DATA			=  2,
		FT_SINO_DATA			=  3,
		FT_MU_DATA				=  4,
		FT_MATRIX_STACK 		=  5,
		FT_IMAGE_STACK			=  6,
		FT_SINO_STACK			=  7,
		FT_CMATRIX_DATA			=  8,
		FT_FLUO_TOMO_IMAGE		=  9,
		FT_ATTENUATION_SINOGRAM = 10,
		FT_DIFFERENCE_SINOGRAM	= 11,
		FT_GENERIC_EDF			= 12,
		FT_GENERIC_VOXEDF		= 13,
		FT_NOT_SPECIFIED		= -1,
		FT_FAMELESS 			= -2,
	};
*/
	
	typedef union
	{
		float			*dt_float;
		double			*dt_double;
		char			*dt_char;
		short			*dt_short;
		int 			*dt_int;
		long			*dt_long;
		unsigned char	*dt_uchar;
		unsigned short	*dt_ushort;
		unsigned int	*dt_uint;
		unsigned long	*dt_ulong;
	} typeconverter_t;
	
	
	typedef std::vector<unsigned long> dims_t;
	typedef std::vector<std::string> text_t;


	//=============================================================== C_HeaderID
	// A class for parsing, storing and writing header IDs.
	//
	class C_HeaderID
	{
		public:
			C_HeaderID();
			C_HeaderID( int curr, int next, int prev );
			C_HeaderID( std::string input );
			
			void init();
			
			void set( int curr, int next, int prev );
			void set_curr( int curr );
			void set_next( int next );
			void set_prev( int prev );
			
			int get_curr() const;
			int get_next() const;
			int get_prev() const;
		
			int str( std::string input );
			std::string str() const;

			bool is_valid() const;
			operator void*() const;
			bool operator!() const;

			std::string print() const;

		private:
			int 			p_curr;
			int 			p_next;
			int 			p_prev;
			std::string 	p_prefix;
			std::string 	p_separator;
			int 			p_fieldsize;
			bool			p_valid;
	};
	std::ostream &operator<<( std::ostream &os, const C_HeaderID &headerID );


	//================================================================ C_Tagpair
	//	A class for managing single tags as used in he edf header.
	//
	class C_Tagpair
	{
		public:
			C_Tagpair();
			C_Tagpair( ns_key::C_Keyword keyword, std::string value, std::string comment="" );
			
			void reset();
			
			void set( ns_key::C_Keyword keyword, std::string value, std::string comment="" );
			void set_keyword( ns_key::C_Keyword keyword );
			void set_value( std::string value );
			void set_value( double value, int precision );
			void set_value( double value );
			void set_value( unsigned long value );
			void set_value( unsigned int value );
			void set_value( int value );
			void set_value( bool value );
			void set_comment( std::string comment );
			
			ns_key::C_Keyword get_keyword() const;
			std::string get_mainKey() const;
			std::string get_value() const;
			std::string get_comment() const;
			
			std::string 	to_string() const;
			double			to_double() const;
			unsigned long	to_ulint() const;
			unsigned int	to_uint() const;
			int 			to_int() const;
			bool			to_bool() const;
			
			std::string   &operator>>( std::string	 &value ) const;
			double		  &operator>>( double		 &value ) const;
			unsigned long &operator>>( unsigned long &value ) const;
			unsigned int  &operator>>( unsigned int  &value ) const;
			int 		  &operator>>( int			 &value ) const;
			bool		  &operator>>( bool 		 &value ) const;
			
			void make_valid();
			bool is_valid() const;
			operator void*() const;
			bool operator!() const;
			
			std::string print() const;

		private:
			ns_key::C_Keyword	p_keyword;	// The part before before assign token.
			std::string p_value;	// The part between assign token and separator.
			std::string p_comment;	// The part between separator and newline.
			bool		p_valid;	// True: No valid has been explicitly set.
	};
	std::ostream &operator<<( std::ostream &os, C_Tagpair tagpair );

	typedef std::vector<C_Tagpair> tags_t;


	//================================================================ C_Taglist
	// A class that manages access to edf header tags.
	//
	class C_Taglist
	{
		public:
			C_Taglist( bool force=true );
			
			void clear();
			void set_force( bool force );
			void set( ns_key::C_Keyword keyword, std::string value );
			void set( ns_key::C_Keyword keyword, const char* value );
			void set( ns_key::C_Keyword keyword, int value );
			void set( ns_key::C_Keyword keyword, unsigned int value );
			void set( ns_key::C_Keyword keyword, unsigned long value );
			void set( ns_key::C_Keyword keyword, double value );
			void set( ns_key::C_Keyword keyword, bool value );
			void set( ns_key::C_Keyword keyword, std::string value, std::string comment );
			void set( ns_key::C_Keyword keyword, const char* value, std::string comment );
			void set( ns_key::C_Keyword keyword, int value, std::string comment );
			void set( ns_key::C_Keyword keyword, unsigned int value, std::string comment );
			void set( ns_key::C_Keyword keyword, unsigned long value, std::string comment );
			void set( ns_key::C_Keyword keyword, double value, std::string comment );
			void set( ns_key::C_Keyword keyword, bool value, std::string comment );
			void set( C_Tagpair tagpair );
			void remove( ns_key::C_Keyword keyword );
			void remove( int index );
			bool get_force() const;
			C_Tagpair get( ns_key::C_Keyword keyword ) const;
			C_Tagpair get( int index ) const;
			int find( C_Tagpair tagpair ) const;
			int find( ns_key::C_Keyword keyword ) const;
			bool is_defined( C_Tagpair tagpair ) const;
			bool is_defined( ns_key::C_Keyword keyword ) const;
			int size() const;
			
			std::vector<std::string> to_vectorlist() const; // acctr<--->Tango
			std::string print() const;

		private:
			tags_t p_tags;
			bool p_force;
	};
	std::ostream &operator<<( std::ostream &os, const C_Taglist &taglist );


	//====================================================================== edf
	// A class for loading and saving edf files.
	//
	class edf 
	{
		public:
			edf( bool verboseval=true );
			edf( std::string name, bool verbose=true );

			void init();

			//static filetype_t convert_str2filetype( std::string filetype_str );
			//static std::string convert_filetype2str( filetype_t filetype );

			//------------------------------------------------------- Set values
			void set_FilenameIn( std::string FilenameIn );
			void set_FilenameOut( std::string FilenameOut );
			void set_HeaderID( const C_HeaderID &headerID );
			void set_Size( unsigned long int size );
			void set_Dim_x( unsigned long dim_x );
			void set_Dim_y( unsigned long dim_y );
			void set_Dim_z( unsigned long dim_z );
			void set_Dim( unsigned long dim_x );
			void set_Dim( unsigned long dim_x, unsigned long dim_y );
			void set_Dim( unsigned long dim_x, unsigned long dim_y, unsigned long dim_z );
			void set_Dimension( int direction, unsigned long dim );
			void set_ByteOrderIn( byteorder_t byteorderIn );
			void set_ByteOrderUser( byteorder_t byteorderOut );
			void set_ByteOrderIn( std::string byteorder_str );
			void set_ByteOrderOutputMode( boom_t ByteOrderOutputMode );
			void set_DataType( datatype_t datatype );
			void set_DataType( std::string datatype_str );
			void set_RealFormat( realformat_t realformat );
			void set_RealFormat( std::string realformat_str );
			void set_ScalingMin( double minval );
			void set_ScalingMax( double maxval );
			void set_FileType( filetype_t filetype );
			//void set_FileType( std::string filetype_str );
			void set_DateOrig( const ns_timer::C_Date &dateOrig );
			void set_DateOrig( std::string dateOrig_str );
			void set_DatePrev( const ns_timer::C_Date &datePrev );
			void set_DatePrev( std::string datePrev_str );
			void set_DateCurr( const ns_timer::C_Date &dateModified );
			void set_DateCurr( std::string dateModified_str );
			void set_FilenameOrig( std::string filenameOrig );
			void set_ImageNumber( int imageNumber );
			void set_VersionNumber( std::string versionNumber );
			void set_ScaleIn( bool scaleIn );
			void set_ScaledFlag( scaled_t scaledFlag );
			void setVerbosity( bool verbose);		//no underscore by matlib convention
			void set_NewFlag( bool newFlag );
			void set_MatlibFlag( bool newFlag );

			//------------------------------------------------------- Get values
			//
			std::string get_FilenameIn() const;
			std::string get_FilenameOut() const;
			C_HeaderID get_HeaderID() const;
			unsigned long get_Size() const;
			unsigned long get_Dim_x() const;
			unsigned long get_Dim_y() const;
			unsigned long get_Dim_z() const;
			unsigned long get_Dimension( int direction ) const;
			int get_NumberOfDimensions() const;
			byteorder_t get_ByteOrderIn() const;
			byteorder_t get_ByteOrderUser() const;
			byteorder_t get_ByteOrderMachine() const;
			byteorder_t get_ByteOrderOut() const;
			boom_t get_ByteOrderOutputMode() const;
			std::string get_ByteOrderIn_str() const;
			std::string get_ByteOrderOut_str() const;
			std::string get_ByteOrder_str( byteorder_t byteorder ) const;
			datatype_t get_DataType() const;
			std::string get_DataType_str() const;
			int get_SizeOfDataType() const;
			bool is_datatype_real() const;
			realformat_t get_RealFormat() const;
			std::string get_RealFormat_str() const;
			double get_ScalingMin() const;
			double get_ScalingMax() const;
			filetype_t get_FileType() const;
			std::string get_FileType_str() const;
			ns_timer::C_Date get_DateOrig() const;
			ns_timer::C_Date get_DatePrev() const;
			ns_timer::C_Date get_DateCurr() const;
			std::string get_DateOrig_str() const;
			std::string get_DatePrev_str() const;
			std::string get_DateCurr_str() const;
			std::string get_FilenameOrig() const;
			std::string get_FilenamePrev() const;
			std::string get_FilenameCurr() const;
			int get_ImageNumber() const;
			std::string get_VersionNumber() const;
			bool get_ScaleIn() const;
			bool get_ScaleOut() const;
			scaled_t get_ScaledFlag() const;
			bool getVerbosity() const;		// no underscore by matlib convention
			bool get_NewFlag() const;
			bool get_MatlibFlag() const;
			
			C_Taglist get_tags_builtin() const;
			C_Taglist get_tags_userdef() const;
			C_Taglist get_tags_all() const;
			bool headerOk( bool readonly=false );
			bool is_complex() const;

			//--------------------------- Functions concerning user defined tags
			//
			void clear_tags();
			void set_tag( ns_key::C_Keyword keyword, std::string value );
			void set_tag( ns_key::C_Keyword keyword, const char* value );
			void set_tag( ns_key::C_Keyword keyword, int value );
			void set_tag( ns_key::C_Keyword keyword, unsigned int value );
			void set_tag( ns_key::C_Keyword keyword, unsigned long value );
			void set_tag( ns_key::C_Keyword keyword, double value );
			void set_tag( ns_key::C_Keyword keyword, bool value );
			void set_tag( ns_key::C_Keyword keyword, std::string value, std::string comment );
			void set_tag( ns_key::C_Keyword keyword, const char* value, std::string comment );
			void set_tag( ns_key::C_Keyword keyword, int value, std::string comment );
			void set_tag( ns_key::C_Keyword keyword, unsigned int value, std::string comment );
			void set_tag( ns_key::C_Keyword keyword, unsigned long value, std::string comment );
			void set_tag( ns_key::C_Keyword keyword, double value, std::string comment );
			void set_tag( ns_key::C_Keyword keyword, bool value, std::string comment );
			void set_tag( C_Tagpair tagpair );
			void remove_tag( ns_key::C_Keyword keyword );
			void remove_tag( int index );
			C_Tagpair get_tag( ns_key::C_Keyword keyword ) const;
			C_Tagpair get_tag( int index ) const;
			int find_tag( C_Tagpair tagpair ) const;
			int find_tag( ns_key::C_Keyword keyword ) const;
			bool is_tag_defined( C_Tagpair tagpair ) const;
			bool is_tag_defined( ns_key::C_Keyword keyword ) const;
			int get_num_tags() const;

			//------------------------------------------- Reading/writing header
			//
			int read_header();
			int read_header( std::string filename );
			int read_header( std::ifstream &file );
			int write_header();
			int write_header( std::string filename );
			int write_header( std::ofstream &file );

			//--------------------------------------- Reading/writing image data
			//
			// Changed interface from 'double *&mat' to 'double *' in all read functions.
			// (JP, 2009-02-10)
			// Changed back to Ô*&Õ because of memory allocation within these functions.
			// (JP, 2009-03-15)
			int read_data( double *&mat, bool flipByteOrder=false );
			int read_data( double *&mat, std::string filename, bool flipByteOrder=false );
			int read_data( double *&mat, std::ifstream &fin, bool flipByteOrder=false );
			int write_data( double *mat, bool flipByteOrder=false );
			int write_data( double *mat, std::string filename, bool flipByteOrder=false );
			int write_data( double *mat, std::ofstream &fout, bool flipByteOrder=false );

			//---------------------------- Reading/writing header and image data
			//
			// Changed interface from 'double *&mat' to 'double *' in all read functions.
			// (JP, 2009-02-10)
			// Changed back to Ô*&Õ because of memory allocation within these functions.
			// (JP, 2009-03-15)
			int read( double *&mat, bool flipByteOrder=false );
			int read( double *&mat, std::string filename, bool flipByteOrder=false );
			int write( double *mat, bool flipByteOrder=false );
			int write( double *mat, std::string filename, bool flipByteOrder=false );

			std::string print() const;

		private:
			std::string 	p_FilenameIn; 	// Image input filename.
			std::string 	p_FilenameOut; 	// Image output filename.
			C_HeaderID		p_HeaderID;		// Format: "EH:%6d:%6d:%6d", current, next, prvious.
			byteorder_t 	p_ByteOrderIn;	// Byte order of binary data in input file.
			byteorder_t 	p_ByteOrderUser;// Output byte order as specified by user.
			boom_t			p_ByteOrderOutputMode;
			datatype_t		p_DataType;		// Data type of data in edf file.
			realformat_t	p_RealFormat;	// Only meaningful for float data types (IEEE, VAX, VMS, CONVEX, etc.).
			dims_t			p_Dims; 		// Number of pixels in directions 1, 2, 3, ...
			unsigned long	p_Size; 		// Number of data bytes (image data only, header not included).
			double			p_ScalingMin;	// Scaling for integer to double conversion.
			double			p_ScalingMax;	// Scaling for integer to double conversion.
			bool			p_scalingMin_defined;
			bool			p_scalingMax_defined;
			filetype_t		p_FileType; 	// matrixdata, sinodata, etc.
			ns_timer::C_Date	p_DateOrig; //
			ns_timer::C_Date	p_DatePrev; //
			ns_timer::C_Date	p_DateCurr; //
			std::string 	p_DateOrig_str; //
			std::string 	p_DatePrev_str; //
			std::string 	p_DateCurr_str; //
			std::string 	p_FilenameOrig;	// Name of the original file when the image was created.
			int 			p_ImageNumber;	// Image number identifying the image.
			std::string 	p_VersionNumber;// Version number of this edf engine.
			std::ios::pos_type	p_BeginData;// Image data begins at this position in current file.
			bool 			p_scaleIn; 		// For reading data from file; depends on edf header.
			scaled_t		p_scaledFlag; 	// For writing data to file; depends on user request.
		                					// and must be treated independently from 'p_scaleIn'.
			bool 			p_verbose;		// Print information to console if true.
			bool			p_newFlag;		// True, if one of tags {TAG_TUD_DATE_ORIG, TAG_TUD_DATE_PREV, TAG_TUD_DATE_CURR} is found. Used to determine, whether byte order of float/double was handled correclty.
			bool			p_matlibFlag;	// True, if file was written by matlib edf engine.

			C_Taglist		p_tags_builtin;		// List of built-in tags (those being interpreted by edf).
			C_Taglist		p_tags_userdef;		// List of user-defined tags (without interpreted ones).
			C_Taglist		p_tags_loaded;		// List of header tags as they have been loaded.
			ns_key::C_AllKeywords	p_allKeywords;	// List of all well-defined keywords, the user is recommended to use.

			int parse_line( std::string input, C_Tagpair &tagpair ) const;
			int interpret_tags();
			unsigned long calculate_headersize( text_t text, int numCharsEndl, int &numSpacePadding ) const;
			unsigned long calculate_datasize() const;
			bool check_consistency( bool readonly=false );

//			void setParams(unsigned long int ny, unsigned long int nx);

			template <class T>		
			int writeTag( text_t &text, ns_key::C_Keyword keyword, T value, std::string comment="" )
			{
				std::ostringstream osst;

				osst << value;
				return writeTagStr( text, keyword, osst.str(), comment );
			};

			template <class T>		
			int writeTagIndex( text_t &text, text_t::size_type index, ns_key::C_Keyword keyword, T value, std::string comment="" )
			{
				std::ostringstream osst;

				osst << value;
				return writeTagStrIndex( text, index, keyword, osst.str(), comment );
			};

			int writeTagStr( text_t &text, ns_key::C_Keyword keyword, std::string value, std::string comment );
			int writeTagStrIndex( text_t &text, text_t::size_type index, ns_key::C_Keyword keyword, std::string value, std::string comment );
	};
	std::ostream &operator<<( std::ostream &os, const edf &file );
}


#endif //edf_h_flag

//==============================================================================
// The following 4 lines are necessary for using meschach's matrix header.
// They are included automatically as soon as <vector> is
// included. This is the reason why there was always trouble
// when including matrixdata.h, depending on which other
// header files had been included before.
// I hope, these 4 lines will solve those problems.
// (by JP on August 28th, 2006)
//
// Removed again, because the reasons for the trouble seem to have vanished
// since edf has been upgraded.
// (by JP on October 12th, 2007)
//
// Note by JF: matrixdata.h does not depend on meschach, anymore.
/*
#ifndef WIN32
	#ifndef _GLIBCXX_VECTOR
	#define _GLIBCXX_VECTOR 1
	    #include <bits/stl_bvector.h>
	#endif
#endif
*/
//==============================================================================
