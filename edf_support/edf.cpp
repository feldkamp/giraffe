//==========================================================================
//	Title: edf.cpp
//
//	Author: Christian G. Schroer
//
//	description: Methods for class ESRF-Data Format
//==========================================================================

#include "edf.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "util.h"

using std::cerr;
using std::cin;
using std::endl;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::ios_base;
using matlibutil::BigEndianByteOrder;
using matlibutil::reverseByteOrder;
using matlibutil::revByteOrder;
using matlibutil::spacepadding;
using matlibutil::toString;
using matlibutil::toInt;
using ns_timer::C_Date;

using namespace ns_key;

namespace ns_edf
{
	//=============================================================== C_HeaderID
	//
	//
	C_HeaderID::C_HeaderID()
	{
		init();
	}
	
	C_HeaderID::C_HeaderID( int curr, int next, int prev )
	{
		set( curr, next, prev );
	}

	C_HeaderID::C_HeaderID( string input )
	{
		str( input );
	}

	void C_HeaderID::init()
	{
		p_curr = 0;
		p_next = 0;
		p_prev = 0;
		p_prefix = HEADER_ID_PREFIX;
		p_separator = HEADER_ID_SEPARATOR;
		p_fieldsize = HEADER_ID_FIELDSIZE;
		p_valid = false;
	}

	void C_HeaderID::set( int curr, int next, int prev )
	{
		set_curr( curr );
		set_next( next );
		set_prev( prev );
	}
	
	void C_HeaderID::set_curr( int curr )
	{
		p_curr = curr;
		p_valid = true;
	}
	
	void C_HeaderID::set_next( int next )
	{
		p_next = next;
		p_valid = true;
	}
	
	void C_HeaderID::set_prev( int prev )
	{
		p_prev = prev;
		p_valid = true;
	}

	int C_HeaderID::get_curr() const
	{
		return p_curr;
	}
	
	int C_HeaderID::get_next() const
	{
		return p_next;
	}
	
	int C_HeaderID::get_prev() const
	{
		return p_prev;
	}

	int C_HeaderID::str( string input )
	{
		string::size_type pos1=0;
		string::size_type pos2=0;
		
		init();
		
		pos1 = input.find( p_separator );
		if ( pos1 == string::npos )
		{
			// Invalid header ID, missing separator.
			return 1;
		}
		
		if ( input.substr( 0, pos1 ) != p_prefix )
		{
			// Missing/invalid suffix.
			return 2;
		}
		
		for ( int i=0; i<3; i++ )
		{
			string value_str;
			istringstream isst;
			
			if ( i<3-1 )
			{
				pos2 = input.find( p_separator, pos1+1 );
				if ( pos2 == string::npos )
				{
					// Missing second separator.
					return 3;
				}
				value_str = input.substr( pos1+1, pos2-pos1-1 );
			}
			else
			{
				value_str = input.substr( pos1+1 );
			}

			if ( (int )value_str.length() != p_fieldsize )
			{
				// Wrong fieldsize.
				return 4;
			}
			isst.str( value_str );
			
			switch ( i )
			{
				case 0: isst >> p_curr; break;
				case 1: isst >> p_next; break;
				case 2: isst >> p_prev; break;
			}
			
			pos1 = pos2;
		}
		
		p_valid = true;
		
		return 0;
	}
	
	string C_HeaderID::str() const
	{
		ostringstream osst;
		
		osst << p_prefix
			 << p_separator << toString( p_curr, p_fieldsize )
			 << p_separator << toString( p_next, p_fieldsize )
			 << p_separator << toString( p_prev, p_fieldsize );
	
		return osst.str();
	}

	bool C_HeaderID::is_valid() const
	{
		return p_valid;
	}
	
	C_HeaderID::operator void*() const
	{
		return is_valid() ? const_cast<C_HeaderID*>(this) : 0;
	}
	
	bool C_HeaderID::operator!() const
	{
		return !is_valid();
	}

	string C_HeaderID::print() const
	{
		ostringstream osst;
		
		osst << "Current header number:  " << get_curr() << endl;
		osst << "Next header number:     " << get_next() << endl;
		osst << "Previous header number: " << get_prev() << endl;
		osst << "Separator: '" << p_separator << "'" << endl;
		osst << "Fieldsize: " << p_fieldsize << endl;
		osst << "Valid: " << ( is_valid() ? "Yes" : "No" ) << endl;
		
		return osst.str();
	}
	
	ostream &operator<<( ostream &os, const C_HeaderID &headerID )
	{
		return os << headerID.print();
	}


	//==========================================================================
	//================================================================ C_Tagpair
	//
	//

	C_Tagpair::C_Tagpair()
	{
		reset();
	}

	C_Tagpair::C_Tagpair( C_Keyword keyword, string value, string comment )
	{
		set( keyword, value, comment );
	}

	void C_Tagpair::reset()
	{
		p_keyword.clear();
		p_value = "";
		p_comment = "";
		p_valid = false;
	}

	void C_Tagpair::set( C_Keyword keyword, string value, string comment )
	{
		set_keyword( keyword );
		set_value( value );
		set_comment( comment );
	}

	void C_Tagpair::set_keyword( C_Keyword keyword )
	{
		p_keyword = keyword;
		make_valid();
	}

	void C_Tagpair::set_value( string value )
	{
		p_value = value;
		make_valid();
	}

	void C_Tagpair::set_value( double value, int precision )
	{
		ostringstream osst;
		
		osst << std::setprecision(precision) << value;
		set_value( osst.str() );
		make_valid();
	}

	void C_Tagpair::set_value( double value )
	{
		ostringstream osst;
		
		osst << value;
		set_value( osst.str() );
		make_valid();
	}

	void C_Tagpair::set_value( unsigned long value )
	{
		ostringstream osst;
		
		osst << value;
		set_value( osst.str() );
	}

	void C_Tagpair::set_value( unsigned int value )
	{
		ostringstream osst;
		
		osst << value;
		set_value( osst.str() );
	}

	void C_Tagpair::set_value( int value )
	{
		ostringstream osst;
		
		osst << value;
		set_value( osst.str() );
	}

	void C_Tagpair::set_value( bool value )
	{
		ostringstream osst;
		
		osst << value;
		set_value( osst.str() );
	}

	void C_Tagpair::set_comment( string comment )
	{
		p_comment = comment;
		make_valid();
	}

	C_Keyword C_Tagpair::get_keyword() const
	{
		return p_keyword;
	}
	
	string C_Tagpair::get_mainKey() const
	{
		return p_keyword.get_main();
	}
	
	string C_Tagpair::get_value() const
	{
		return p_value;
	}
	
	string C_Tagpair::get_comment() const
	{
		return p_comment;
	}

	// Exactly the same as C_Tagpair::get_value().
	// For symmetry resons only...
	string C_Tagpair::to_string() const
	{
		return p_value;
	}
	
	double C_Tagpair::to_double() const
	{
		double value;
		
		*this >> value;
		return value;
	}
	
	unsigned long C_Tagpair::to_ulint() const
	{
		unsigned long value;
		
		*this >> value;
		return value;
	}
	
	unsigned int C_Tagpair::to_uint() const
	{
		unsigned int value;
		
		*this >> value;
		return value;
	}
	
	int C_Tagpair::to_int() const
	{
		int value;
		
		*this >> value;
		return value;
	}
	
	bool C_Tagpair::to_bool() const
	{
		bool value;
		
		*this >> value;
		return value;
	}
	
	string &C_Tagpair::operator>>( string &value ) const
	{
		return value = get_value();
	}
	
	double &C_Tagpair::operator>>( double &value ) const
	{
		istringstream isst;
		
		isst.str( get_value() );
		isst >> value;
		
		return value;
	}
	
	unsigned long &C_Tagpair::operator>>( unsigned long &value ) const
	{
		istringstream isst;
		
		isst.str( get_value() );
		isst >> value;
		
		return value;
	}
	
	unsigned int &C_Tagpair::operator>>( unsigned int &value ) const
	{
		istringstream isst;
		
		isst.str( get_value() );
		isst >> value;
		
		return value;
	}
	
	int &C_Tagpair::operator>>( int &value ) const
	{
		istringstream isst;
		
		isst.str( get_value() );
		isst >> value;
		
		return value;
	}
	
	bool &C_Tagpair::operator>>( bool &value ) const
	{
		istringstream isst;
		
		isst.str( get_value() );
		isst >> value;
		
		return value;
	}
	
	void C_Tagpair::make_valid()
	{
		p_valid = true;
	}
	
	bool C_Tagpair::is_valid() const
	{
		return p_valid;
	}
	
	C_Tagpair::operator void*() const
	{
		return is_valid() ? const_cast<C_Tagpair*>(this) : 0;
	}
	
	bool C_Tagpair::operator!() const
	{
		return !is_valid();
	}

	string C_Tagpair::print() const
	{
		ostringstream osst;
		
		osst << "Keywords: " << get_keyword() << endl;
		osst << "Value:    " << get_value() << endl;
		osst << "Comment:  " << get_comment() << endl;
		osst << "Valid:    " << ( is_valid() ? "yes" : "no" ) << endl;
		
		return osst.str();
	}
	
	ostream &operator<<( ostream &os, C_Tagpair tagpair )
	{
		return os << tagpair.print();
	}
	

	//==========================================================================
	//================================================================ C_Taglist
	//
	//
	C_Taglist::C_Taglist( bool force )
	{
		p_force = force;
	}

	//-------------------------------------------------------------------- clear
	//
	//
	void C_Taglist::clear()
	{
		p_tags.clear();
	}

	//---------------------------------------------------------------- set_force
	//
	//
	void C_Taglist::set_force( bool force )
	{
		p_force = force;
	}
	
	//---------------------------------------------------------------------- set
	//
	//
	void C_Taglist::set( C_Keyword keyword, string value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, const char* value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( string(value) );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, int value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, unsigned int value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, unsigned long value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, double value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, bool value )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, string value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, const char* value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( string(value) );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, int value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, unsigned int value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, unsigned long value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, double value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Keyword keyword, bool value, string comment )
	{
		C_Tagpair tagpair;
		tagpair.set_keyword( keyword );
		tagpair.set_value( value );
		tagpair.set_comment( comment );
		set( tagpair );
	}
	
	void C_Taglist::set( C_Tagpair tagpair )
	{
		// If someone thinks that performance is crucial, he/she may improve
		// performance by checking existance via 'index==-1' instead of the
		// defined() function. But I think no taglist will be large enough,
		// so checking first and finding then is quite ok and more safe.
		// (Safe, because returning data type of find may be changed to
		// unsigned, etc.)

		if ( !is_defined( tagpair ) )
		{
			// Keyword does not exist in taglist. Add new tagpair.
			p_tags.push_back( tagpair );
		}
		else if ( get_force() )
		{
			int index;
			
			// Keyword exists and force is set true. Replace old tagpair.
			index = find( tagpair );
			p_tags[index].set_value( tagpair.get_value() );
			p_tags[index].set_comment( tagpair.get_comment() );
		}
		else
		{
			// Keyword exists and force is toggled off. Do nothing.
		}
	}

	//------------------------------------------------------------------- remove
	//
	//
	void C_Taglist::remove( C_Keyword keyword )
	{
		if ( is_defined( keyword ) )
		{
			p_tags.erase( p_tags.begin() + find( keyword ) );
		}
	}
	
	//------------------------------------------------------------------- remove
	//
	//
	void C_Taglist::remove( int index )
	{
		if ( index >= 0 && index < size() )
		{
			p_tags.erase( p_tags.begin() + index );
		}
	}

	//---------------------------------------------------------------- get_force
	//
	//
	bool C_Taglist::get_force() const
	{
		return p_force;
	}
	
	//---------------------------------------------------------------------- get
	//
	//
	C_Tagpair C_Taglist::get( C_Keyword keyword ) const
	{
		C_Tagpair tagpair;

		// Concerning performance see comment in C_Taglist::set().
		if ( is_defined( keyword ) )
		{
			tagpair = p_tags[find(keyword)];
		}
		
		return tagpair;
	}
	
	//------------------------------------------------------------------ get_tag
	//
	//
	C_Tagpair C_Taglist::get( int index ) const
	{
		C_Tagpair tagpair;
		
		if ( index >= 0 && index < size() )
		{
			tagpair = p_tags[index];
		}
		
		return tagpair;
	}

	//--------------------------------------------------------------------- find
	//
	//
	int C_Taglist::find( C_Tagpair tagpair ) const
	{
		return find( tagpair.get_keyword() );
	}

	//--------------------------------------------------------------------- find
	//
	//
	int C_Taglist::find( C_Keyword keyword ) const
	{
		int index;

		index = -1;
		for ( int i=0; i<size(); i++ )
		{
			if ( get(i).get_keyword() == keyword )
			{
				index = i;
				break;
			}
		}

		return index;
	}
	
	//------------------------------------------------------------------ defined
	//
	//
	bool C_Taglist::is_defined( C_Tagpair tagpair ) const
	{
		return is_defined( tagpair.get_keyword() );
	}
	
	//------------------------------------------------------------------ defined
	//
	//
	bool C_Taglist::is_defined( C_Keyword keyword ) const
	{
		bool defined;

		defined = false;
		for ( int i=0; i<size(); i++ )
		{
			if ( get(i).get_keyword() == keyword )
			{
				defined = true;
				break;
			}
		}

		return defined;
	}

	//--------------------------------------------------------------------- size
	//
	//
	int C_Taglist::size() const
	{
		return (int )p_tags.size();
	}

	//-------------------------------------------------------------------- print
	//
	// This function converts a taglist to a vector of strings.
	// The format is <key1 value1 comment1 key2 value2 comment2 ...>.
	// key1, key2, ... are the main keys; the key aliases will be skipped.
	// The intention is to make it possible to pass the Taglist through
	// Tango (Corba) interfaces as can be seen in counter.cpp, mca.cpp and
	// camera.cpp.
	//
	vector<string> C_Taglist::to_vectorlist() const
	{
		vector<string> vectorlist;

		for ( int i=0; i<size(); i++ )
		{
			vectorlist.push_back( get(i).get_mainKey() );
			vectorlist.push_back( get(i).get_value() );
			vectorlist.push_back( get(i).get_comment() );
		}
	
		return vectorlist;
	}
	
	//-------------------------------------------------------------------- print
	//
	//
	string C_Taglist::print() const
	{
		ostringstream osst;
		
		osst << "Force: " << ( get_force() ? "Yes" : "No" ) << endl;
		osst << "Defined tags: " << endl;
		osst << "{" << endl;
		for ( tags_t::const_iterator i = p_tags.begin(); i != p_tags.end(); i++ )
		{
			osst << *i;
		}
		osst << "}" << endl;
		
		return osst.str();
	}
	
	//---------------------------------------------------------------- operator<<
	//
	//
	ostream &operator<<( ostream &os, const C_Taglist &taglist )
	{
		return os << taglist.print();
	}


	//==========================================================================
	//====================================================================== edf
	//
	//
	edf::edf( bool verbose )
	{
		init();
		setVerbosity( verbose );	// verbose output to console
	}

	//---------------------------------------------------------------------- edf
	//
	//
	edf::edf( string filename, bool verbose )
	{
		init();
		set_FilenameIn( filename );
		read_header();
		setVerbosity( verbose );	// verbose output to console
	}

	//--------------------------------------------------------------------- init
	//
	//
	void edf::init()
	{
		p_FilenameIn			= "";
		p_FilenameOut			= "";
		p_HeaderID.init();
		p_ByteOrderIn 			= get_ByteOrderMachine();
		p_ByteOrderUser			= get_ByteOrderMachine();
		p_ByteOrderOutputMode	= BOOM_MACHINE;
		p_DataType				= DT_UNSIGNED_SHORT;
		p_RealFormat 			= REAL_FORMAT_IEEE;
		p_Dims.clear();
		p_Size					= 0;
		p_ScalingMin			= 0;
		p_ScalingMax			= 0;
		p_scalingMin_defined	= false;
		p_scalingMax_defined	= false;
		p_FileType				= FTS_GENERIC_EDF;
		p_DateOrig.init();
		p_DatePrev.init();
		p_DateCurr.init();
		p_DateOrig_str			= "";
		p_DatePrev_str			= "";
		p_DateCurr_str			= "";
		p_FilenameOrig			= "";
		p_ImageNumber			= 1;
		p_VersionNumber 		= EDF_VERSION_NUMBER;
		p_BeginData 			= 0;
		p_scaleIn				= false;
		p_scaledFlag			= SF_SCALED;
		p_verbose				= true;
		p_newFlag				= false;
		p_matlibFlag			= false;
		p_tags_builtin.clear();
		p_tags_userdef.clear();
		p_tags_loaded.clear();
	}

	//---------------------------------------------------------------------- set
	//
	//
	void edf::set_FilenameIn( string FilenameIn )
	{
		p_FilenameIn = FilenameIn;
	}

	void edf::set_FilenameOut( string FilenameOut )
	{
		p_FilenameOut = FilenameOut;
	}

	void edf::set_HeaderID( const C_HeaderID &HeaderID )
	{
		p_HeaderID = HeaderID;
	}

	void edf::set_Size(unsigned long int Size)
	{
		p_Size = Size;
	}

	void edf::set_Dim_x( unsigned long dim_x )
	{
		set_Dimension( 0, dim_x );
	}

	void edf::set_Dim_y( unsigned long dim_y )
	{
		set_Dimension( 1, dim_y );
	}

	void edf::set_Dim_z( unsigned long dim_z )
	{
		set_Dimension( 2, dim_z );
	}

	void edf::set_Dim( unsigned long dim_x )
	{
		set_Dim_x( dim_x );
	}
	
	void edf::set_Dim( unsigned long dim_x, unsigned long dim_y )
	{
		set_Dim_x( dim_x );
		set_Dim_y( dim_y );
	}
	
	void edf::set_Dim( unsigned long dim_x, unsigned long dim_y, unsigned long dim_z )
	{
		set_Dim_x( dim_x );
		set_Dim_y( dim_y );
		set_Dim_z( dim_z );
	}
	
	void edf::set_Dimension( int direction, unsigned long dim )
	{
		while ( (int )p_Dims.size() <= direction )
		{
			p_Dims.push_back(1);
		}
		p_Dims[direction] = dim;
	}

	void edf::set_ByteOrderIn( byteorder_t ByteOrderIn )
	{
		p_ByteOrderIn = ByteOrderIn;
	}

	void edf::set_ByteOrderUser( byteorder_t ByteOrderUser )
	{
		p_ByteOrderUser = ByteOrderUser;
	}

	void edf::set_ByteOrderIn( string byteorder_str )
	{
		if ( byteorder_str == HIGH_BYTE_FIRST_STR )
		{
			p_ByteOrderIn = BO_HIGH_BYTE_FIRST;
		}
		else if ( byteorder_str == LOW_BYTE_FIRST_STR )
		{
			p_ByteOrderIn = BO_LOW_BYTE_FIRST;
		}
		else if ( byteorder_str == NOT_SPECIFIED_STR_TUD ||
				  byteorder_str == NOT_SPECIFIED_STR_ESRF )
		{
			p_ByteOrderIn = BO_NOT_SPECIFIED;
		}
		else
		{
			p_ByteOrderIn = BO_FAMELESS;
		}
	}

	void edf::set_ByteOrderOutputMode( boom_t ByteOrderOutputMode )
	{
		p_ByteOrderOutputMode = ByteOrderOutputMode;
	}

	void edf::set_DataType( datatype_t DataType )
	{
		p_DataType = DataType;
	}

	void edf::set_DataType( string datatype_str )
	{
		if ( datatype_str == FLOAT_SINGLE_STR ||
			 datatype_str == OLD_FLOAT_STR )
		{
			p_DataType = DT_FLOAT;
		}
		else if ( datatype_str == FLOAT_DOUBLE_STR ||
				  datatype_str == OLD_DOUBLE_STR )
		{
			p_DataType = DT_DOUBLE;
		}
		else if ( datatype_str == BYTE_STR )
		{
			p_DataType = DT_BYTE;
		}
		else if ( datatype_str == SHORT_STR ||
				  datatype_str == OLD_SHORT_STR )
		{
			p_DataType = DT_SHORT;
		}
		else if ( datatype_str == INTEGER_STR )
		{
			p_DataType = DT_INTEGER;
		}
		else if ( datatype_str == LONG_STR ||
				  datatype_str == OLD_LONG_STR )
		{
			p_DataType = DT_LONG;
		}
		else if ( datatype_str == UNSIGNED_BYTE_STR )
		{
			p_DataType = DT_UNSIGNED_BYTE;
		}
		else if ( datatype_str == UNSIGNED_SHORT_STR ||
				  datatype_str == OLD_UNSIGNED_SHORT_STR )
		{
			p_DataType = DT_UNSIGNED_SHORT;
		}
		else if ( datatype_str == UNSIGNED_INTEGER_STR )
		{
			p_DataType = DT_UNSIGNED_INTEGER;
		}
		else if ( datatype_str == UNSIGNED_LONG_STR ||
				  datatype_str == OLD_UNSIGNED_LONG_STR )
		{
			p_DataType = DT_UNSIGNED_LONG;
		}
		else if ( datatype_str == COMPLEX_SINGLE_STR )
		{
			p_DataType = DT_COMPLEX_FLOAT;
		}
		else if ( datatype_str == COMPLEX_DOUBLE_STR )
		{
			p_DataType = DT_COMPLEX_DOUBLE;
		}
		else if ( datatype_str == COMPLEX_BYTE_STR )
		{
			p_DataType = DT_COMPLEX_BYTE;
		}
		else if ( datatype_str == COMPLEX_SHORT_STR )
		{
			p_DataType = DT_COMPLEX_SHORT;
		}
		else if ( datatype_str == COMPLEX_INTEGER_STR )
		{
			p_DataType = DT_COMPLEX_INTEGER;
		}
		else if ( datatype_str == COMPLEX_LONG_STR )
		{
			p_DataType = DT_COMPLEX_LONG;
		}
		else if ( datatype_str == NOT_SPECIFIED_STR_TUD ||
				  datatype_str == NOT_SPECIFIED_STR_ESRF )
		{
			p_DataType = DT_NOT_SPECIFIED;
		}
		else
		{
			p_DataType = DT_FAMELESS;
		}
	}
		
	void edf::set_RealFormat( realformat_t RealFormat )
	{
		p_RealFormat = RealFormat;
	}

	void edf::set_RealFormat( string realformat_str )
	{
		if ( realformat_str == REAL_FORMAT_IEEE_STR )
		{
			p_RealFormat = REAL_FORMAT_IEEE;
		}
		else if ( realformat_str == REAL_FORMAT_VAX_STR )
		{
			p_RealFormat = REAL_FORMAT_VAX;
		}
		else if ( realformat_str == REAL_FORMAT_VMS_STR )
		{
			p_RealFormat = REAL_FORMAT_VMS;
		}
		else if ( realformat_str == REAL_FORMAT_CONVEX_STR )
		{
			p_RealFormat = REAL_FORMAT_CONVEX;
		}
		else if ( realformat_str == NOT_SPECIFIED_STR_TUD ||
				  realformat_str == NOT_SPECIFIED_STR_ESRF )
		{
			p_RealFormat = REAL_FORMAT_NOT_SPECIFIED;
		}
		else
		{
			p_RealFormat = REAL_FORMAT_FAMELESS;
		}
	}

	void edf::set_ScalingMin( double ScalingMin )
	{
		p_ScalingMin = ScalingMin;
		p_scalingMin_defined = true;
	}

	void edf::set_ScalingMax( double ScalingMax )
	{
		p_ScalingMax = ScalingMax;
		p_scalingMax_defined = true;
	}

	void edf::set_FileType( filetype_t FileType )
	{
		p_FileType = FileType;
	}
/*
	void edf::set_FileType( string filetype_str )
	{
		if ( filetype_str == FT_VECTOR_DATA_STR )
		{
			p_FileType = FT_VECTOR_DATA;
		}
		else if ( filetype_str == FT_MATRIX_DATA_STR )
		{
			p_FileType = FT_MATRIX_DATA;
		}
		else if ( filetype_str == FT_IMAGE_DATA_STR )
		{
			p_FileType = FT_IMAGE_DATA;
		}
		else if ( filetype_str == FT_SINO_DATA_STR )
		{
			p_FileType = FT_SINO_DATA;
		}
		else if ( filetype_str == FT_MU_DATA_STR )
		{
			p_FileType = FT_MU_DATA;
		}
		else if ( filetype_str == FT_MATRIX_STACK_STR )
		{
			p_FileType = FT_MATRIX_STACK;
		}
		else if ( filetype_str == FT_IMAGE_STACK_STR )
		{
			p_FileType = FT_IMAGE_STACK;
		}
		else if ( filetype_str == FT_SINO_STACK_STR )
		{
			p_FileType = FT_SINO_STACK;
		}
		else if ( filetype_str == FT_CMATRIX_DATA_STR )
		{
			p_FileType = FT_CMATRIX_DATA;
		}
		else if ( filetype_str == FT_FLUO_TOMO_IMAGE_STR )
		{
			p_FileType = FT_FLUO_TOMO_IMAGE;
		}
		else if ( filetype_str == FT_ATTENUATION_SINOGRAM_STR )
		{
			p_FileType = FT_ATTENUATION_SINOGRAM;
		}
		else if ( filetype_str == FT_DIFFERENCE_SINOGRAM_STR )
		{
			p_FileType = FT_DIFFERENCE_SINOGRAM;
		}
		else if ( filetype_str == FT_GENERIC_EDF_STR )
		{
			p_FileType = FT_GENERIC_EDF;
		}
		else if ( filetype_str == FT_GENERIC_VOXEDF_STR )
		{
			p_FileType = FT_GENERIC_VOXEDF;
		}
		else if ( filetype_str == NOT_SPECIFIED_STR_TUD ||
				  filetype_str == NOT_SPECIFIED_STR_ESRF )
		{
			p_FileType = FT_NOT_SPECIFIED;
		}
		else
		{
			p_FileType = FT_FAMELESS;
		}
	}
*/
	void edf::set_DateOrig( const C_Date &DateOrig )
	{
		p_DateOrig = DateOrig;
	}
	
	void edf::set_DateOrig( string DateOrig_str )
	{
		p_DateOrig_str = DateOrig_str;
	}

	void edf::set_DatePrev( const C_Date &DatePrev )
	{
		p_DatePrev = DatePrev;
	}

	void edf::set_DatePrev( string DatePrev_str )
	{
		p_DatePrev_str = DatePrev_str;
	}

	void edf::set_DateCurr( const C_Date &DateCurr )
	{
		p_DateCurr = DateCurr;
	}

	void edf::set_DateCurr( string DateCurr_str )
	{
		p_DateCurr_str = DateCurr_str;
	}

	void edf::set_FilenameOrig( string FilenameOrig )
	{
		p_FilenameOrig = FilenameOrig;
	}

	void edf::set_ImageNumber( int ImageNumber )
	{
		p_ImageNumber = ImageNumber;
	}

	void edf::set_VersionNumber( string VersionNumber )
	{
		p_VersionNumber = VersionNumber;
	}

	void edf::set_ScaleIn(bool scaleIn)
	{
		p_scaleIn = scaleIn;
	}

	void edf::set_ScaledFlag( scaled_t scaledFlag)
	{
		p_scaledFlag = scaledFlag;
	}

	void edf::setVerbosity(bool verbose)
	{
		p_verbose = verbose;
	}

	void edf::set_NewFlag(bool newFlag)
	{
		p_newFlag = newFlag;
	}

	void edf::set_MatlibFlag(bool matlibFlag)
	{
		p_matlibFlag = matlibFlag;
	}
	
	//---------------------------------------------------------------------- get
	//
	//
	string edf::get_FilenameIn() const
	{
		return p_FilenameIn;
	}

	string edf::get_FilenameOut() const
	{
		return p_FilenameOut;
	}

	C_HeaderID edf::get_HeaderID() const
	{
		return p_HeaderID;
	}

	unsigned long int edf::get_Size() const
	{
		return p_Size;
	}

	unsigned long edf::get_Dim_x() const
	{
		return get_Dimension(0);
	}

	unsigned long edf::get_Dim_y() const
	{
		return get_Dimension(1);
	}

	unsigned long edf::get_Dim_z() const
	{
		return get_Dimension(2);
	}

	unsigned long edf::get_Dimension( int direction ) const
	{ 
		unsigned long result;
		
		if ( direction >= (int )p_Dims.size() )
		{
			result = 0;
		}
		else
		{
			result = p_Dims[direction];
		}
		return result;
	}

	int edf::get_NumberOfDimensions() const
	{
		return (int )p_Dims.size();
	}

	byteorder_t edf::get_ByteOrderIn() const
	{
		return p_ByteOrderIn;
	}

	byteorder_t edf::get_ByteOrderMachine() const
	{
		return BigEndianByteOrder() ? BO_HIGH_BYTE_FIRST : BO_LOW_BYTE_FIRST;
	}

	byteorder_t edf::get_ByteOrderUser() const
	{
		return p_ByteOrderUser;
	}
	
	byteorder_t edf::get_ByteOrderOut() const
	{
		byteorder_t ByteOrderOut;
		
		switch ( get_ByteOrderOutputMode() )
		{
			case BOOM_INPUT:	ByteOrderOut = get_ByteOrderIn();	   break;
			case BOOM_MACHINE:	ByteOrderOut = get_ByteOrderMachine(); break;
			case BOOM_USER: 	ByteOrderOut = get_ByteOrderUser();    break;
			default:
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::get_ByteOrderOut(). ";
				cerr <<	"Invalid byteorder output mode '" << get_ByteOrderOutputMode() << "'." << endl;
				ByteOrderOut = get_ByteOrderMachine();
		}
		return ByteOrderOut;
	}

	boom_t edf::get_ByteOrderOutputMode() const
	{
		return p_ByteOrderOutputMode;
	}

	string edf::get_ByteOrderIn_str() const
	{
		return get_ByteOrder_str( get_ByteOrderIn() );
	}
	
	string edf::get_ByteOrderOut_str() const
	{
		return get_ByteOrder_str( get_ByteOrderOut() );
	}
	
	string edf::get_ByteOrder_str( byteorder_t byteorder ) const
	{
		string byteorder_str;

		switch ( byteorder )
		{
			case BO_HIGH_BYTE_FIRST:	byteorder_str = HIGH_BYTE_FIRST_STR;	break;
			case BO_LOW_BYTE_FIRST:		byteorder_str = LOW_BYTE_FIRST_STR;		break;
			case BO_NOT_SPECIFIED:		byteorder_str = NOT_SPECIFIED_STR_ESRF; break;
			case BO_FAMELESS:			byteorder_str = FAMELESS_STR;			break;
			default:
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::get_ByteOrder_str(). ";
				cerr <<	"Invalid byteorder '" << byteorder << "'." << endl;
				byteorder_str = FAMELESS_STR;
		}

		return byteorder_str;
	}

	datatype_t edf::get_DataType() const
	{
		return p_DataType;
	}

	string edf::get_DataType_str() const
	{
		string datatype_str;

		switch ( get_DataType() )
		{
			case DT_FLOAT:				datatype_str = FLOAT_SINGLE_STR;		break;
			case DT_DOUBLE: 			datatype_str = FLOAT_DOUBLE_STR;		break;

			case DT_BYTE:				datatype_str = BYTE_STR; 				break;
			case DT_SHORT:				datatype_str = SHORT_STR;				break;
			case DT_INTEGER:			datatype_str = INTEGER_STR;				break;
			case DT_LONG:				datatype_str = LONG_STR; 				break;

			case DT_UNSIGNED_BYTE:		datatype_str = UNSIGNED_BYTE_STR;		break;
			case DT_UNSIGNED_SHORT:		datatype_str = UNSIGNED_SHORT_STR;		break;
			case DT_UNSIGNED_INTEGER:	datatype_str = UNSIGNED_INTEGER_STR;	break;
			case DT_UNSIGNED_LONG:		datatype_str = UNSIGNED_LONG_STR;		break;

			case DT_COMPLEX_FLOAT:		datatype_str = COMPLEX_SINGLE_STR;		break;
			case DT_COMPLEX_DOUBLE:		datatype_str = COMPLEX_DOUBLE_STR;		break;
			case DT_COMPLEX_BYTE:		datatype_str = COMPLEX_BYTE_STR;		break;
			case DT_COMPLEX_SHORT:		datatype_str = COMPLEX_SHORT_STR;		break;
			case DT_COMPLEX_INTEGER:	datatype_str = COMPLEX_INTEGER_STR;		break;
			case DT_COMPLEX_LONG:		datatype_str = COMPLEX_LONG_STR;		break;

			case DT_NOT_SPECIFIED:		datatype_str = NOT_SPECIFIED_STR_ESRF;	break;
			case DT_FAMELESS:			datatype_str = FAMELESS_STR;			break;
			default:
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::get_DataType_str(). ";
				cerr <<	"Invalid data type '" << get_DataType() << "'." << endl;
				datatype_str = FAMELESS_STR;
		}

		return datatype_str;
	}

	//------------------------------------------------------- get_SizeOfDataType
	//
	// Returns the storage size (in bytes) of given data types in the edf file.
	// In contrast to previous implementations of edf, the data size is not
	// retrieved from the c language sizeof operator, but it is hardcoded.
	// This is, because we are not interested in the storage in main memory
	// but in the storage of the file on disk. Thus we need a value independent
	// from the c language and standardized for edf.
	//
	int edf::get_SizeOfDataType() const
	{
		int size;

		switch (get_DataType())
		{
			case DT_FLOAT:				size =  4; break;
			case DT_DOUBLE: 			size =  8; break;

			case DT_BYTE:				size =  1; break;
			case DT_SHORT:				size =  2; break;
			case DT_INTEGER:			size =  4; break;
			case DT_LONG:				size =  4; break;

			case DT_UNSIGNED_BYTE:		size =  1; break;
			case DT_UNSIGNED_SHORT:		size =  2; break;
			case DT_UNSIGNED_INTEGER:	size =  4; break;
			case DT_UNSIGNED_LONG:		size =  4; break;

			case DT_COMPLEX_FLOAT:		size =  8; break;
			case DT_COMPLEX_DOUBLE:		size = 16; break;
			case DT_COMPLEX_BYTE:		size =  2; break;
			case DT_COMPLEX_SHORT:		size =  4; break;
			case DT_COMPLEX_INTEGER:	size =  8; break;
			case DT_COMPLEX_LONG:		size =  8; break;

			case DT_NOT_SPECIFIED:		size =  0; break;
			case DT_FAMELESS:			size =  0; break;

			default:
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::getSizeOfDataType(). ";
				cerr <<	"Invalid data type '" << get_DataType() << "'." << endl;
				size = 0;
		}

		return size;
	}

	bool edf::is_datatype_real() const
	{
		bool is_real;
	
		switch (get_DataType())
		{
			case DT_FLOAT:				is_real = true;  break;
			case DT_DOUBLE: 			is_real = true;  break;

			case DT_BYTE:				is_real = false; break;
			case DT_SHORT:				is_real = false; break;
			case DT_INTEGER:			is_real = false; break;
			case DT_LONG:				is_real = false; break;

			case DT_UNSIGNED_BYTE:		is_real = false; break;
			case DT_UNSIGNED_SHORT:		is_real = false; break;
			case DT_UNSIGNED_INTEGER:	is_real = false; break;
			case DT_UNSIGNED_LONG:		is_real = false; break;

			case DT_COMPLEX_FLOAT:		is_real = true;  break;
			case DT_COMPLEX_DOUBLE:		is_real = true;  break;
			case DT_COMPLEX_BYTE:		is_real = false; break;
			case DT_COMPLEX_SHORT:		is_real = false; break;
			case DT_COMPLEX_INTEGER:	is_real = false; break;
			case DT_COMPLEX_LONG:		is_real = false; break;

			case DT_NOT_SPECIFIED:		is_real = false; break;
			case DT_FAMELESS:			is_real = false; break;

			default:
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::is_datatype_real(). ";
				cerr <<	"Invalid data type '" << get_DataType() << "'." << endl;
				is_real = false;
		}
		
		return is_real;
	}
	
	realformat_t edf::get_RealFormat() const
	{
		return p_RealFormat;
	}

	string edf::get_RealFormat_str() const
	{
		string realformat_str;
		
		switch (get_RealFormat())
		{
			case REAL_FORMAT_IEEE:			realformat_str = REAL_FORMAT_IEEE_STR;	 break;
			case REAL_FORMAT_VAX:			realformat_str = REAL_FORMAT_VAX_STR;	 break;
			case REAL_FORMAT_VMS:			realformat_str = REAL_FORMAT_VMS_STR;	 break;
			case REAL_FORMAT_CONVEX:		realformat_str = REAL_FORMAT_CONVEX_STR; break;
			case REAL_FORMAT_NOT_SPECIFIED: realformat_str = NOT_SPECIFIED_STR_ESRF; break;
			case REAL_FORMAT_FAMELESS: 		realformat_str = FAMELESS_STR;			 break;
			default:
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::get_RealFormat_str(). ";
				cerr <<	"Invalid real format '" << get_RealFormat() << "'." << endl;
				realformat_str = FAMELESS_STR;
		}
	
		return realformat_str;
	}

	double edf::get_ScalingMin() const
	{
		return p_ScalingMin;
	}

	double edf::get_ScalingMax() const
	{
		return p_ScalingMax;
	}

	filetype_t edf::get_FileType() const
	{
		return p_FileType;
	}

	string edf::get_FileType_str() const
	{
		return p_FileType;
	}

/*
	filetype_t edf::convert_str2filetype( string filetype_str )
	{

		if ( filetype_str == FTS_VECTOR_DATA )
		{
			return FT_VECTOR_DATA;
		}
		else if ( filetype_str == FTS_MATRIX_DATA )
		{
			return FT_MATRIX_DATA;
		}
		else if ( filetype_str == FTS_IMAGE_DATA )
		{
			return FT_IMAGE_DATA;
		}
		else if ( filetype_str == FTS_SINO_DATA )
		{
			return FT_SINO_DATA;
		}
		else if ( filetype_str == FTS_MU_DATA )
		{
			return FT_MU_DATA;
		}
		else if ( filetype_str == FTS_MATRIX_STACK )
		{
			return FT_MATRIX_STACK;
		}
		else if ( filetype_str == FTS_IMAGE_STACK )
		{
			return FT_IMAGE_STACK;
		}
		else if ( filetype_str == FTS_SINO_STACK )
		{
			return FT_SINO_STACK;
		}
		else if ( filetype_str == FTS_CMATRIX_DATA )
		{
			return FT_CMATRIX_DATA;
		}
		else if ( filetype_str == FTS_DIFFERENCE_SINOGRAM )
		{
			return FT_DIFFERENCE_SINOGRAM;
		}
		else if ( filetype_str == FTS_GENERIC_EDF )
		{
			return FT_GENERIC_EDF;
		}
		else if ( filetype_str == FTS_GENERIC_VOXEDF )
		{
			return FT_GENERIC_VOXEDF;
		}
		else if ( filetype_str == FTS_FLUO_TOMO_IMAGE )
		{
			return FT_FLUO_TOMO_IMAGE;
		}
		else if ( filetype_str == FTS_ATTENUATION_SINOGRAM )
		{
			return FT_ATTENUATION_SINOGRAM;
		}
		else if ( filetype_str == NOT_SPECIFIED_STR_ESRF )
		{
			return FT_NOT_SPECIFIED;
		}
		else
		{
			return FT_FAMELESS;
		}
	}

	string edf::convert_filetype2str( filetype_t filetype )
	{
		string fileType_str;
		
		switch ( filetype )
		{
			case FT_VECTOR_DATA 		: fileType_str = FTS_VECTOR_DATA		  ; break;
			case FT_MATRIX_DATA 		: fileType_str = FTS_MATRIX_DATA		  ; break;
			case FT_IMAGE_DATA			: fileType_str = FTS_IMAGE_DATA			  ; break;
			case FT_SINO_DATA			: fileType_str = FTS_SINO_DATA			  ; break;
			case FT_MU_DATA 			: fileType_str = FTS_MU_DATA			  ; break;
			case FT_MATRIX_STACK		: fileType_str = FTS_MATRIX_STACK		  ; break;
			case FT_IMAGE_STACK 		: fileType_str = FTS_IMAGE_STACK		  ; break;
			case FT_SINO_STACK			: fileType_str = FTS_SINO_STACK			  ; break;
			case FT_CMATRIX_DATA		: fileType_str = FTS_CMATRIX_DATA		  ; break;
			case FT_FLUO_TOMO_IMAGE 	: fileType_str = FTS_DIFFERENCE_SINOGRAM  ; break;
			case FT_ATTENUATION_SINOGRAM: fileType_str = FTS_GENERIC_EDF 		  ; break;
			case FT_DIFFERENCE_SINOGRAM : fileType_str = FTS_GENERIC_VOXEDF		  ; break;
			case FT_GENERIC_EDF 		: fileType_str = FTS_FLUO_TOMO_IMAGE 	  ; break;
			case FT_GENERIC_VOXEDF		: fileType_str = FTS_ATTENUATION_SINOGRAM ; break;
			case FT_NOT_SPECIFIED		: fileType_str = NOT_SPECIFIED_STR_ESRF   ; break; 
			case FT_FAMELESS			: fileType_str = FAMELESS_STR;			  ; break;
			default:
				// The user should never see this message.
				cerr << "ERROR. In matrixstack::convert_filetype2string() - ";
				cerr <<	"Invalid filetype '" << filetype << "'." << endl;
				fileType_str = FAMELESS_STR;
		}
		
		return fileType_str;
	}
*/
	C_Date edf::get_DateOrig() const
	{
		return p_DateOrig;
	}
	
	C_Date edf::get_DatePrev() const
	{
		return p_DatePrev;
	}

	C_Date edf::get_DateCurr() const
	{
		return p_DateCurr;
	}

	string edf::get_DateOrig_str() const
	{
		string date;
		
		if ( get_FilenameIn() == "" )
		{
			date = C_Date().get_date_now().print();
		}
		else
		{
			if ( p_DateOrig_str == "" )
			{
				// Please do not remove the incommented lines!
				
				date = C_Date().get_date_now().print();
				//date = "";
				//date = NOT_SPECIFIED_STR_ESRF;
				//date = FAMELESS_STR;
			}
			else
			{
				date = p_DateOrig_str;
			}
		}
		
		return date;
	}
	
	string edf::get_DatePrev_str() const
	{
		string date;
		
		if ( get_FilenameIn() == "" )
		{
			date = C_Date().get_date_now().print();
		}
		else
		{
			if ( p_DatePrev_str == "" )
			{
				// Please do not remove the incommented lines!
				
				date = C_Date().get_date_now().print();
				//date = "";
				//date = NOT_SPECIFIED_STR_ESRF;
				//date = FAMELESS_STR;
			}
			else
			{
				date = p_DatePrev_str;
			}
		}
		
		return date;
	}

	string edf::get_DateCurr_str() const
	{
		return C_Date().get_date_now().print();
	}

	string edf::get_FilenameOrig() const
	{
		string filename;
		
		if ( get_FilenameIn() == "" )
		{
			filename = get_FilenameOut();
		}
		else
		{
			if ( p_FilenameOrig == "" )
			{
				// Please do not remove the incommented lines!
				
				filename = get_FilenameIn();
				//filename = "";
				//filename = NOT_SPECIFIED_STR_ESRF;
				//filename = FAMELESS_STR;
			}
			else
			{
				filename = p_FilenameOrig;
			}
		}
				
		return filename;
	}

	string edf::get_FilenamePrev() const
	{
		string filename;
		
		if ( get_FilenameIn() == "" )
		{
			filename = get_FilenameOut();
		}
		else
		{
			filename = get_FilenameIn();
		}
				
		return filename;
	}

	string edf::get_FilenameCurr() const
	{
		return get_FilenameOut();
	}

	int edf::get_ImageNumber() const
	{
		return p_ImageNumber;
	}

	string edf::get_VersionNumber() const
	{
		return p_VersionNumber;
	}

	bool edf::get_ScaleIn() const
	{
		return p_scaleIn;
	}

	bool edf::get_ScaleOut() const
	{
		bool scaleOut=0;
		
		switch ( get_ScaledFlag() )
		{
			case SF_UNSCALED:	scaleOut = false;			break;
			case SF_SCALED: 	scaleOut = true;			break;
			case SF_RETAIN: 	scaleOut = get_ScaleIn();	break;
		}
		
		return scaleOut;
	}

	scaled_t edf::get_ScaledFlag() const
	{
		return p_scaledFlag;
	}

	bool edf::getVerbosity() const
	{
		return p_verbose;
	}

	bool edf::get_NewFlag() const
	{
		return p_newFlag;
	}

	bool edf::get_MatlibFlag() const
	{
		return p_matlibFlag;
	}

	C_Taglist edf::get_tags_builtin() const
	{
		return p_tags_builtin;
	}
	
	C_Taglist edf::get_tags_userdef() const
	{
		return p_tags_userdef;
	}
	
	C_Taglist edf::get_tags_all() const
	{
		return p_tags_loaded;
	}

	bool edf::headerOk( bool readonly )
	{
		return check_consistency( readonly );
	}

	bool edf::is_complex() const
	{
		if ( ( get_DataType() == DT_COMPLEX_FLOAT   ) ||
			 ( get_DataType() == DT_COMPLEX_DOUBLE  ) ||
			 ( get_DataType() == DT_COMPLEX_BYTE    ) ||
			 ( get_DataType() == DT_COMPLEX_SHORT   ) ||
			 ( get_DataType() == DT_COMPLEX_INTEGER ) ||
			 ( get_DataType() == DT_COMPLEX_LONG    ) )
		{
			return true;
		}
		return false;
	}
	
	//-------------------------------------------------------------- clear_tags
	//
	//
	void edf::clear_tags()
	{
		p_tags_userdef.clear();
	}
			
	void edf::set_tag( C_Keyword keyword, string value )
	{
		//cout << "A string: value=" << value << endl;
		p_tags_userdef.set( keyword, value, "" );
	}
	
	// This overloaded function is necessary, because literal string constants
	// would be interpreted as bool values insteead as strings.
	void edf::set_tag( C_Keyword keyword, const char* value )
	{
		//cout << "A const char*: \nkeyword=" << keyword << ", value=" << value << endl;
		p_tags_userdef.set( keyword, string(value), "" );
	}
	
	void edf::set_tag( C_Keyword keyword, int value )
	{
		//cout << "A int: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), "" );
	}
	
	void edf::set_tag( C_Keyword keyword, unsigned int value )
	{
		//cout << "A unsigned int: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), "" );
	}
	
	void edf::set_tag( C_Keyword keyword, unsigned long value )
	{
		//cout << "A unsigned long: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), "" );
	}
	
	void edf::set_tag( C_Keyword keyword, double value )
	{
		//cout << "A double: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), "" );
	}
	
	void edf::set_tag( C_Keyword keyword, bool value )
	{
		//cout << "A bool: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), "" );
	}
	
	void edf::set_tag( C_Keyword keyword, string value, string comment )
	{
		//cout << "B string: value=" << value << endl;
		p_tags_userdef.set( keyword, value, comment );
	}
	
	// This overloaded function is necessary, because literal string constants
	// would be interpreted as bool values insteead as strings.
	void edf::set_tag( C_Keyword keyword, const char* value, string comment )
	{
		//cout << "B const char*: value=" << value << endl;
		p_tags_userdef.set( keyword, string(value), comment );
	}
	
	void edf::set_tag( C_Keyword keyword, int value, string comment )
	{
		//cout << "B int: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), comment );
	}
	
	void edf::set_tag( C_Keyword keyword, unsigned int value, string comment )
	{
		//cout << "B unsigned int: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), comment );
	}
	
	void edf::set_tag( C_Keyword keyword, unsigned long value, string comment )
	{
		//cout << "B unsigned long: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), comment );
	}
	
	void edf::set_tag( C_Keyword keyword, double value, string comment )
	{
		//cout << "B double: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), comment );
	}
	
	void edf::set_tag( C_Keyword keyword, bool value, string comment )
	{
		//cout << "B bool: value=" << value << endl;
		ostringstream osst;
		osst << value;
		p_tags_userdef.set( keyword, osst.str(), comment );
	}
	
	void edf::set_tag( C_Tagpair tagpair )
	{
		p_tags_userdef.set( tagpair );
	}
	
	void edf::remove_tag( C_Keyword keyword )
	{
		p_tags_userdef.remove( keyword );
	}
	
	void edf::remove_tag( int index )
	{
		p_tags_userdef.remove( index );
	}
	
	C_Tagpair edf::get_tag( C_Keyword keyword ) const
	{
		return p_tags_userdef.get( keyword );
	}
	
	C_Tagpair edf::get_tag( int index ) const
	{
		return p_tags_userdef.get( index );
	}
	
	int edf::find_tag( C_Tagpair tagpair ) const
	{
		return p_tags_userdef.find( tagpair );
	}
	
	int edf::find_tag( C_Keyword keyword ) const
	{
		return p_tags_userdef.find( keyword );
	}
	
	bool edf::is_tag_defined( C_Tagpair tagpair ) const
	{
		return p_tags_userdef.is_defined( tagpair );
	}
	
	bool edf::is_tag_defined( C_Keyword keyword ) const
	{
		return p_tags_userdef.is_defined( keyword );
	}
	
	int edf::get_num_tags() const
	{
		return p_tags_userdef.size();
	}
	


	//-------------------------------------------------------------- read_header
	//
	//
	int edf::read_header()
	{
		return read_header( get_FilenameIn() );
	}
	
	//-------------------------------------------------------------- read_header
	//
	//
	int edf::read_header( string filename )
	{
		ifstream fin;
		int retval;
		
		set_FilenameIn( filename );

		fin.open( filename.c_str() );
		if ( !fin )
		{
			cout << "ERROR. Could not open file '" << filename << "' for reading header." << endl;
			return 1;
		}
		retval = read_header( fin );
		fin.close();
		
		return retval;
	}
	
	//-------------------------------------------------------------- read_header
	//
	//
	int edf::read_header( ifstream &fin )
	{
		string line;
 		bool found;
		int count;
		int retval;

		p_tags_loaded.clear();
		p_tags_builtin.clear();
		p_tags_userdef.clear();

		if ( !fin )
		{
			return 1;
		}

		found = false;
		while ( getline( fin, line ) )
		{
			if ( line.find( HEADER_TOKEN_BEGIN ) != string::npos )
			{
				found = true;
				break;
			}
			else if ( line != "" )
			{
				cout << "ERROR. This file seems not to be in edf format." << endl;
				cout << " Between the begin of the file and the header token '";
				cout << HEADER_TOKEN_BEGIN << "' must be only endline characters." << endl;
				return 2;
			}
		}
		
		if ( !found )
		{
			cout << "ERROR. Missing header token '" << HEADER_TOKEN_BEGIN << "'." << endl;;
			return 3;
		}

		found = false;
		count = 1;
		while ( getline( fin, line ) && !found )
		{
			C_Tagpair tagpair;
			int retval;
	
			count++;
			retval = parse_line( line, tagpair );

			switch (retval)
			{
				case 0:
				{
					if ( tagpair )
					{
						p_tags_loaded.set( tagpair );
					}
					else
					{
						// Invalid tagpair (resulting from empty lines or lines
						// with only whitespaces or separators) will be ignored.
					}
				} break;
			
				case 1: // The end of the header has been reached.
				{
					std::ios::pos_type filepos;
					
					filepos = fin.tellg();
					p_BeginData = filepos;
					if ( filepos % HEADER_BYTES_MODULO != 0 )
					{
						cout << "WARNING. Header size is not an integer multiple ";
						cout << "of " << HEADER_BYTES_MODULO << "." << endl;
					}
					found = true;
				} break;
				
				case 2:
				{
					cout << "WARNING. In line " << count << ". Missing separator '";
					cout << HEADER_TOKEN_SEPARATOR << "'. Skip this line." << endl;
				} break;

				case 3:
				{
					cout << "WARNING. In line " << count << ". Missing assign token '";
					cout << HEADER_TOKEN_ASSIGN << "'. Skip this line." << endl;
				} break;

				case 4:
				{
					cout << "WARNING. In line " << count << ". Missing keyword before '";
					cout << HEADER_TOKEN_ASSIGN << "'. Skip this line." << endl;
				} break;

				default:
				{
					cout << "WARNING. In line " << count << ". Unknown error ";
					cout << retval << ". Skip this line." << endl;
				} break;
			}
		}

		if ( !found )
		{
			return 5;
		}
		
		//cout << "1> dimensions: " << get_NumberOfDimensions() << endl;
		retval = interpret_tags();
		//cout << "2> dimensions: " << get_NumberOfDimensions() << endl;
		if ( retval )
		{
			return 6;
		}

		return 0;
	}
	
	//------------------------------------------------------------- write_header
	//
	//
	int edf::write_header()
	{
		return write_header( get_FilenameOut() );
	}
	
	//------------------------------------------------------------- write_header
	//
	//
	int edf::write_header( string filename )
	{
		ofstream fout;

		set_FilenameOut( filename );
		
		fout.open( filename.c_str(), ios_base::binary );
		if ( !fout )
		{
			cout << "ERROR. Could not open file '" << filename << "' for writing header." << endl;
			return 1;
		}
		//cout << "INFO: File '" << filename << "' was opened for writing header." << endl;
		write_header( fout );
		//cout << "INFO: Header has been written." << endl;
		fout.close();
		
		return 0;
	}

	//------------------------------------------------------------- write_header
	//
	//
	int edf::write_header( ofstream &fout )
	{
		text_t text;
		unsigned int headersize;
		int numCharsEndl=0;
		int numSpacePadding=0;
		
		if ( !fout )
		{
			return 1;
		}
		//cout << "DEBUGINFO: edf::write_header() - Enter." << endl;

		// Write the header in datastructure text_t in order to have random access.
		// Random access is needed to fill in the headersize value after it will have
		// been calculated. But for calculating headersize, we need to write the
		// header first... (See function calculate_headersize().)
		text.clear();
		text.push_back( HEADER_TOKEN_BEGIN );
		writeTag( text, TAG_HEADER_ID, get_HeaderID().str() );
		writeTag( text, TAG_TUD_HEADER_SIZE, "" );
		writeTag( text, TAG_TUD_HEADER_LINES, "" );
		writeTag( text, TAG_SIZE, calculate_datasize() );
		for ( int i=0; i<get_NumberOfDimensions(); i++ )
		{
			C_Keyword keyword( TAG_DIM.get_main() + toString(i+1) );
			writeTag( text, keyword, get_Dimension(i) );
		}
		writeTag( text, TAG_BYTE_ORDER, get_ByteOrderOut_str() );
		writeTag( text, TAG_DATA_TYPE,	get_DataType_str() );
		if ( is_datatype_real() )
		{
			writeTag( text, TAG_REAL_FORMAT, get_RealFormat_str() );
		}
		if ( get_ScaleOut() )
		{
			writeTag( text, TAG_TUD_SCALING_MIN, get_ScalingMin() );
			writeTag( text, TAG_TUD_SCALING_MAX, get_ScalingMax() );
		}
		writeTag( text, TAG_TUD_FILE_TYPE, get_FileType_str() );
		writeTag( text, TAG_TUD_DATE_ORIG, get_DateOrig_str() );
		writeTag( text, TAG_TUD_DATE_PREV, get_DatePrev_str() );
		writeTag( text, TAG_TUD_DATE_CURR, get_DateCurr_str() );
		writeTag( text, TAG_TUD_FILENAME_ORIG, get_FilenameOrig() );
		writeTag( text, TAG_TUD_FILENAME_PREV, get_FilenamePrev() );
		writeTag( text, TAG_TUD_FILENAME_CURR, get_FilenameCurr() );
		writeTag( text, TAG_IMAGE, get_ImageNumber() );
		writeTag( text, TAG_VERSION_NUMBER, EDF_VERSION_NUMBER );
		//cout << "DEBUGINFO: edf::write_header() - AAA" << endl;
		// Add user defined header keywords.
		for ( int i=0; i<get_num_tags(); i++ )
		{
			C_Tagpair tagpair;

			tagpair = get_tag(i);
			writeTag( text, tagpair.get_keyword(), tagpair.get_value(), tagpair.get_comment() );
		}

		text.push_back( HEADER_TOKEN_END );

		//cout << "DEBUGINFO: edf::write_header() - BBB" << endl;

		// Fill in the number of lines of the header.
		writeTagIndex( text, 3, TAG_TUD_HEADER_LINES, (unsigned int )text.size() );

		// Now we are able to calculate the headersize.
		numCharsEndl = (int )string( HEADER_TOKEN_ENDLINE ).length();
		headersize = (unsigned int)calculate_headersize( text, numCharsEndl, numSpacePadding );

		//cout << "DEBUGINFO: edf::write_header() - numSpacePadding = " << numSpacePadding << endl;

		// Fill in the headersize value.
		writeTagIndex( text, 2, TAG_TUD_HEADER_SIZE, headersize );

		// Fill in as many spaces as needed for a header size that is an interger
		// multiple of 2^10. Do not forget the closing curled bracket '}'.
		text[text.size()-1] = string( numSpacePadding, ' ' ) + HEADER_TOKEN_END;

		// Write the header ascii text to the given file..
		for ( text_t::size_type i=0; i<text.size(); i++ )
		{
			fout << text[i] << HEADER_TOKEN_ENDLINE;
		}
		
		return 0;
	}

	//---------------------------------------------------------------- read_data
	//
	//
	int edf::read_data( double *&mat, bool flipByteOrder )
	{
		return read_data( mat, get_FilenameIn(), flipByteOrder );
	}

	//---------------------------------------------------------------- read_data
	//
	//
	int edf::read_data( double *&mat, string filename, bool flipByteOrder )
	{
		ifstream fin;
		int retval;
		
		set_FilenameIn( filename );

		fin.open( filename.c_str(), ios_base::binary );
		if ( !fin )
		{
			cout << "ERROR. Could not open file '" << filename << "' for reading data." << endl;
			return 1;
		}
		fin.seekg( p_BeginData );
		retval = read_data( mat, fin, flipByteOrder );
		fin.close();

		p_scalingMin_defined = false;
		p_scalingMax_defined = false;

		return retval;
	}

	//---------------------------------------------------------------- read_data
	//
	//
	int edf::read_data( double *&mat, ifstream &fin, bool flipByteOrder )
	{
		string complex_str;
		unsigned long nz=0;
		unsigned long ny=0;
		unsigned long nx=0;
		unsigned long numBytes;
		unsigned long numPixels1;
		unsigned long numPixels2;
		int sizeDataType1;	// For complex types: size of real or imag part
		int sizeDataType2;	// For complex types: size of the whole number
		double value;
		bool reverse_ByteOrder;
		int retval=0;
		
		if ( !fin )
		{
			return 1;
		}

		nx = ( ( get_Dim_x() > 0 ) ? get_Dim_x() : 1 );
		ny = ( ( get_Dim_y() > 0 ) ? get_Dim_y() : 1 );
		nz = ( ( get_Dim_z() > 0 ) ? get_Dim_z() : 1 );

		set_Dim_x( nx ); 
		set_Dim_y( ny );
		set_Dim_z( nz );

		numPixels1 = nx*ny*nz;
		numPixels2 = ( is_complex() ? 2*nx*ny*nz : nx*ny*nz );
		numBytes = calculate_datasize();
		sizeDataType2 = get_SizeOfDataType();
		sizeDataType1 = ( is_complex() ? sizeDataType2/2 : sizeDataType2 );
		complex_str = ( is_complex() ? " complex " : " " );

		if ( mat )
		{
			delete []mat;
			mat=0;
		}
		
		mat = new double[numPixels2];

		if ( ( get_ByteOrderIn() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) ||
			 ( get_ByteOrderIn() == BO_LOW_BYTE_FIRST  ) && (  BigEndianByteOrder() ) )
		{
			reverse_ByteOrder = !flipByteOrder;
		}
		else
		{
			reverse_ByteOrder = flipByteOrder;
		}

		// The following section has been introduced for backward compability.
		// In the past, there was an error when writing float or double data.
		// The byteorder specified in the edf header was the opposite of the
		// byteorder the data was written in the data section of the edf file.
		// This error has been eliminated, but we should still be able to read
		// old edf files with errornous byteorder specification.
		// The funtion get_NewFlag returns true if one of the folloeng headaer
		// tags have been found when reading the header:
		// TAG_TUD_DATE_ORIG, TAG_TUD_DATE_PREV, TAG_TUD_DATE_CURR.
		// (2009-09-02, JP)
		// Pay attention to the fact that complex data types need
		// not to be considered, bacause they have been introduced after the
		// bug removal.
		// (2011-05-19, JP)
		// Bad bug: The absence of the tags TAG_TUD_DATE_ORIG, etc. does NOT
		// necessarily mean, that the edf was written by an old tomo edf engine.
		// The file could have also been written by some other software (like
		// PyMca, for example ;-)). In this case, the byte order must not be
		// changed. To distinguish between the two cases, I introduced the
		// flag 'matlibFlag'.
		if ( ( get_DataType() == DT_FLOAT ||
			   get_DataType() == DT_DOUBLE ) &&
			   !get_NewFlag() && get_MatlibFlag() )
		{
			reverse_ByteOrder = !reverse_ByteOrder;
			cout << "WARNING. You are loading an edf file which had been written" << endl;
			cout << "with an errornous edf engine (tomo version 9.0 or older)." << endl;
		}
		else
		{
			//cout << "INFO: NewFlag has been set." << endl;
		}

		switch ( get_DataType() )
		{
			case DT_FLOAT:
			case DT_COMPLEX_FLOAT:
			{
				float *data;
				
				data = new float[numPixels2];
				
				fin.read( (char *)data, numBytes );
				if ( !fin )
				{
					cout << "WARNING. In edf::read_data() - Could not read "
						<< numBytes << " from file." << endl;
				}

				if ( getVerbosity() )
				{
					cout << numPixels1 << complex_str << "floats read" << endl;
				}

				if ( reverse_ByteOrder )
				{
					revByteOrder( data, numPixels2 );
				}

				for (unsigned long i = 0; i < numPixels2; i++)
				{
					mat[i] = data[i];
				}
				delete []data;		
			}
			break;

			case DT_DOUBLE:
			case DT_COMPLEX_DOUBLE:
			{
				double *data; // For aesthetic reasons...
				
				data = mat;
				
				fin.read( (char *)data, numBytes );
				if ( !fin )
				{
					cout << "WARNING. In edf::read_data() - Could not read "
						<< numBytes << " from file." << endl;
				}
				
				if ( getVerbosity() )
				{
					cout << numPixels1 << complex_str << "doubles read" << endl;
				}
				
				if ( reverse_ByteOrder )
				{
					revByteOrder( data, numPixels2 );
				}
			}
			break;

			case DT_BYTE:				// fallthru
			case DT_SHORT:				// fallthru
			case DT_INTEGER:			// fallthru
			case DT_LONG:				// fallthru
			case DT_COMPLEX_BYTE:		// fallthru
			case DT_COMPLEX_SHORT:		// fallthru
			case DT_COMPLEX_INTEGER:	// fallthru
			case DT_COMPLEX_LONG:		// fallthru
			{
				char *data;
				double range;
				double offset;
				double shift;

				data = new char[ sizeDataType1 * numPixels2 ];
				fin.read( (char *)data, numBytes );
				if ( !fin )
				{
					cout << "WARNING. In edf::read_data() - Could not read "
						<< numBytes << " from file." << endl;
				}
				
				if ( getVerbosity() )
				{
					cout << numPixels1 << " (signed) " << sizeDataType2
						<< "-byte" << complex_str << "integers read" << endl;
				}

				if ( flipByteOrder )
				{
					// Attention: reversing the byte order of complex data types
					// must NOT flip the order of real part and imaginary part.
					reverseByteOrder( data, sizeDataType1, (int)numPixels2 );
				}

				if ( get_ScaleIn() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
					shift  = std::pow( 256.0, sizeDataType1 ) / 2.0;
				}
				else
				{
					range  = 1;
					offset = 0;
					shift  = 0;
				}

				for (unsigned long i = 0; i < numPixels2; i++)
				{
					value = 0;
					for ( int k = 0; k < sizeDataType1; k++ )
					{
						unsigned char byte;

						if ( get_ByteOrderIn() == BO_LOW_BYTE_FIRST  )
						{
							byte = (unsigned char )data[ i*sizeDataType1 + k ];
						}
						else
						{
							byte = (unsigned char )data[ i*sizeDataType1 + sizeDataType1 - k - 1 ];
						}
						
						value += byte * std::pow( 256.0, k );
					}
					if ( value >= std::pow( 256.0, sizeDataType1 ) / 2 )
					{
						value -= std::pow( 256.0, sizeDataType1 );
					}
					mat[i] = ( value + shift ) * range + offset;
				}
				delete []data;
			}
			break;

			case DT_UNSIGNED_BYTE:		// fallthru
			case DT_UNSIGNED_SHORT:		// fallthru
			case DT_UNSIGNED_INTEGER:	// fallthru
			case DT_UNSIGNED_LONG:		// fallthru
			{
				char *data;
				double range;
				double offset;

				data = new char[ sizeDataType1 * numPixels2 ];
				fin.read( data, numBytes );
				if ( !fin )
				{
					cout << "WARNING. In edf::read_data() - Could not read "
						<< numBytes << " from file." << endl;
				}

				if ( getVerbosity() )
				{
					cout << numPixels1 << " (unsigned) " << sizeDataType2
						<< "-byte" << complex_str << "integers read" << endl;
				}

				if ( flipByteOrder )
				{
					reverseByteOrder( data, sizeDataType1, (int)numPixels2 );
				}

				if ( get_ScaleIn() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
				}
				else
				{
					range  = 1;
					offset = 0;
				}

				for (unsigned long i = 0; i < numPixels2; i++)
				{
					value = 0;
					for ( int k = 0; k < sizeDataType1; k++ )
					{
						unsigned char byte;

						if ( get_ByteOrderIn() == BO_LOW_BYTE_FIRST  )
						{
							byte = (unsigned char )data[ i*sizeDataType1 + k ];
						}
						else
						{
							byte = (unsigned char )data[ i*sizeDataType1 + sizeDataType1 - k - 1 ];
						}
						
						value += byte * std::pow( 256.0, k );
					}
					mat[i] = value * range + offset;
				}
				delete []data;
			}
			break;

			case DT_NOT_SPECIFIED:
			{

			}
			break;

			case DT_FAMELESS:
			{

			}
			break;

			default:
			{
				// The user should never see this message.
				cout << "ERROR. In edf::get_datatype_str() - "
					<< "Invalid data type '" << get_DataType() << "'." << endl;
				set_DataType( DT_NOT_SPECIFIED );
				retval = 1;
			}
			break;
		}

		return retval;
	}

	//--------------------------------------------------------------- write_data
	//
	//
	int edf::write_data( double *mat, bool flipByteOrder )
	{
		return write_data( mat, get_FilenameOut(), flipByteOrder );
	}

	//--------------------------------------------------------------- write_data
	//
	//
	int edf::write_data( double *mat, string filename, bool flipByteOrder )
	{
		ofstream fout;
		int retval;
		
		set_FilenameOut( filename );

		fout.open( filename.c_str(), ios_base::app | ios_base::binary );
		if ( !fout )
		{
			cout << "ERROR. Could not open file '" << filename << "' for writing data." << endl;
			return 1;
		}
		retval = write_data( mat, fout, flipByteOrder );
		fout.close();
		
		return retval;
	}

	//--------------------------------------------------------------- write_data
	//
	//
	int edf::write_data( double *mat, ofstream &fout, bool flipByteOrder )
	{
		string complex_str;
		int retval;
		int sizeDataType1;
		int sizeDataType2;
		bool reverse_ByteOrder;
		unsigned long numBytes;
		unsigned long numPixels1;
		unsigned long numPixels2;
		double epsilon; // Added to double values before they are converted to integers
						// in order to ensure the scaled maximum is 2^n-1 instead of 2^n-2.

		if ( get_ScaleOut() )
		{
			if ( !p_scalingMin_defined || !p_scalingMax_defined )
			{
				cout << "ERROR. Scaled flag is set, but min/max are not ";
				cout << "defined. No image data written." << endl;
				return 2;
			}
		}
		
		numBytes = calculate_datasize();
		numPixels1 = get_Dim_y()*get_Dim_x();
		numPixels2 = ( is_complex() ? 2*numPixels1 : numPixels1 );
		
		sizeDataType2 = get_SizeOfDataType();
		sizeDataType1 = ( is_complex() ? sizeDataType2/2 : sizeDataType2 );

		complex_str = ( is_complex() ? " complex " : " " );

		epsilon = 0.1;
		retval = 0;

		if ( ( get_ByteOrderOut() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) ||
			 ( get_ByteOrderOut() == BO_LOW_BYTE_FIRST  ) && (  BigEndianByteOrder() ) )
		{
			reverse_ByteOrder = !flipByteOrder;
		}
		else
		{
			reverse_ByteOrder = flipByteOrder;
		}
		/*
		cout << "reverse ByteOrder: " << (reverse_ByteOrder ? "yes" : "no") << endl;
		cout << "Scale out? "<< (get_ScaleOut() ? "Yes" : "No") << endl;
		*/
		switch ( get_DataType() )
		{	
			case DT_FLOAT:
			case DT_COMPLEX_FLOAT:
			{
				float *data;
				data = new float[numPixels2];

				for (unsigned long i = 0; i < numPixels2; i++)
				{
						data[i] = (float )mat[i];
				}
				/*
				cout << endl;
				cout << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << endl;
				cout << data[4] << " " << data[5] << " " << data[6] << " " << data[7] << endl;
				cout << endl;
				*/
				if ( reverse_ByteOrder )
				{
					// revByteOrder() is substituted by reverseByteOrder().
					// If there should be problems, change it back (JP, 2009-09-02)
					revByteOrder( data, numPixels2 );
					//reverseByteOrder( data, sizeDataType1, numPixels2 );
				}
				/*
				cout << endl;
				cout << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << endl;
				cout << data[4] << " " << data[5] << " " << data[6] << " " << data[7] << endl;
				cout << endl;
				*/

				fout.write( (char *)data, numBytes );
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					cout << numPixels1 << complex_str << "floats written to " << get_FilenameOut() << endl;
				}

				delete []data;	
			}
			break;

			case DT_DOUBLE:
			case DT_COMPLEX_DOUBLE:
			{
				double *data; // For aesthetic reasons...

				// Performance considerations.
				// Two possibilities:
				// a) mat --> revByteOrder(mat) --> write(mat) --> revByteOrder(mat)
				// b) mat --> copy(data, mat) --> revByteOrder(data) --> write(data)
				// In my opinion b) would be faster, but I did no tests until now.
				// So a) stays implemented.
				// Update: b) would be a bad idea for very large arrays, because
				// we could get out of memory...
				// (JP Oct 2007)
				
				data = mat;
				
				if ( reverse_ByteOrder )
				{
					// revByteOrder() is substituted by reverseByteOrder().
					// If there should be problems, change it back (JP, 2009-09-02)
					revByteOrder( data, numPixels2 );
					//reverseByteOrder( data, sizeDataType1, numPixels2 );
				}

				fout.write( (char *)data, numBytes );
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					cout << numPixels1 << complex_str << "doubles written to " << get_FilenameOut() << endl;
				}
				
				if ( reverse_ByteOrder )
				{
					// revByteOrder() is substituted by reverseByteOrder().
					// If there should be problems, change it back (JP, 2009-09-02)
					revByteOrder( data, numPixels2 );
					//reverseByteOrder( data, sizeDataType1, numPixels2 );
				}
			}
			break;

			case DT_BYTE:
			case DT_COMPLEX_BYTE:
			{
				char *data;
				double range;
				double offset;
				double shift;

				data = new char[numPixels2];

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
					shift  = std::pow( 256.0, sizeDataType1 ) / 2.0;
				}
				else
				{
					range  = 1;
					offset = 0;
					shift  = 0;
				}
				
				/*
				cout << "get_ScalingMax()=" << get_ScalingMax() << endl;
				cout << "get_ScalingMin()=" << get_ScalingMin() << endl;
				cout << "range=" << range << endl;
				cout << "offset=" << offset << endl;
				cout << "shift=" << shift << endl;
				cout << endl;
				cout << (int )mat[0] << " " << (int )mat[1] << " " << (int )mat[2] << " " << (int )mat[3] << endl;
				cout << (int )mat[4] << " " << (int )mat[5] << " " << (int )mat[6] << " " << (int )mat[7] << endl;
				cout << endl;
				*/
				
				for (unsigned long i = 0; i < numPixels2; i++)
				{
					if ( range != 0 )
					{
						data[i] = (char )floor( (mat[i]-offset)/range-shift+epsilon );
					}
					else
					{
						data[i] = 0;
					}
				}
				/*	
				cout << endl;
				cout << (int )data[0] << " " << (int )data[1] << " " << (int )data[2] << " " << (int )data[3] << endl;
				cout << (int )data[4] << " " << (int )data[5] << " " << (int )data[6] << " " << (int )data[7] << endl;
				cout << endl;
				*/

				// No byte order for this data type...

				fout.write( (char *)data, numBytes );
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numPixels1 << complex_str << "bytes written to "
						<< get_FilenameOut() << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_SHORT:
			case DT_COMPLEX_SHORT:
			{
				short *data;
				double range;
				double offset;
				double shift;

				data = new short[numPixels2];

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
					shift  = std::pow( 256.0, sizeDataType1 ) / 2.0;
				}
				else
				{
					range  = 1;
					offset = 0;
					shift  = 0;
				}

				for (unsigned long i = 0; i < numPixels2; i++)
				{
					if ( range != 0 )
					{
						data[i] = (short )floor( (mat[i]-offset)/range-shift+epsilon );
					}
					else
					{
						data[i] = 0;
					}
				}		

				if ( reverse_ByteOrder )
				{
					revByteOrder( data, numPixels2 );
				}

				fout.write( (char *)data, numBytes );
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numPixels1 << complex_str << "shorts written to "
						<< get_FilenameOut() << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_INTEGER:			// fallthru
			case DT_LONG:				// fallthru
			case DT_COMPLEX_INTEGER:	// fallthru
			case DT_COMPLEX_LONG:		// fallthru
			{
				int *data_int;
				long *data_long;
				double range;
				double offset;
				double shift;

				if ( sizeof(int) == 4 )
				{
					data_int = new int[numPixels2];
				}
				else if ( sizeof(long) == 4 )
				{
					data_long = new long[numPixels2];
				}
				else
				{
					cout << "ERROR. In edf::write_data() - Neither 'int' nor 'long int' are of size 4 byte." << endl;
					return 3;
				}

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
					shift  = std::pow( 256.0, sizeDataType1 ) / 2.0;
				}
				else
				{
					range  = 1;
					offset = 0;
					shift  = 0;
				}
				for (unsigned long i = 0; i < numPixels2; i++)
				{
					if ( range != 0 )
					{
						if ( sizeof(int) == 4 )
						{
							data_int[i] = (int )floor( (mat[i]-offset)/range-shift+epsilon);
						}
						else
						{
							data_long[i] = (long )floor( (mat[i]-offset)/range-shift+epsilon);
						}
					}
					else
					{
						if ( sizeof(int) == 4 )
						{
							data_int[i] = 0;
						}
						else
						{
							data_long[i] = 0;
						}
					}
				}		

				if ( reverse_ByteOrder )
				{
					if ( sizeof(int) == 4 )
					{
						revByteOrder( data_int, numPixels2 );
					}
					else
					{
						revByteOrder( data_long, numPixels2 );
					}
				}

				if ( sizeof(int) == 4 )
				{
					fout.write( (char *)data_int, numBytes );
				}
				else
				{
					fout.write( (char *)data_long, numBytes );
				}
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numPixels1 << complex_str << "longs written to "
						<< get_FilenameOut() << " (" << un << "scaled)" << endl;
				}

				if ( sizeof(int) == 4 )
				{
					delete []data_int;
				}
				else
				{
					delete []data_long;
				}
			}
			break;
			
			case DT_UNSIGNED_BYTE:
			{
				unsigned char *data;
				double range;
				double offset;

				data = new unsigned char[numPixels2];

				if ( get_ScaleOut() )
				{
					range = ( get_ScalingMax() - get_ScalingMin() )/( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
				}
				else
				{
					range = 1;
					offset = 0;
				}

				for (unsigned long i = 0; i < numPixels2; i++)
				{
					if ( range != 0 )
					{
						data[i] = (unsigned char )floor( (mat[i]-offset)/range+epsilon );
					}
					else
					{
						data[i] = 0;
					}
				}

				// No byte order for this data type...

				fout.write( (char *)data, numBytes );
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numPixels1 << complex_str << "unsigned bytes written to "
						<< get_FilenameOut() << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_UNSIGNED_SHORT:
			{
				unsigned short *data;
				double range;
				double offset;

				data = new unsigned short[numPixels2];

				if ( get_ScaleOut() )
				{
					range = ( get_ScalingMax() - get_ScalingMin() )/( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
				}
				else
				{
					range = 1;
					offset = 0;
				}

				for (unsigned long i = 0; i < numPixels2; i++)
				{
					if ( range != 0 )
					{
						data[i] = (unsigned short )floor( (mat[i]-offset)/range+epsilon );
					}
					else
					{
						data[i] = 0;
					}
				}

				if ( reverse_ByteOrder )
				{
					revByteOrder( data, numPixels2 );
				}

				fout.write( (char *)data, numBytes );
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numPixels1 << complex_str << "unsigned shorts written to "
						<< get_FilenameOut() << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_UNSIGNED_INTEGER:	// fallthru
			case DT_UNSIGNED_LONG:		// fallthru
			{
				unsigned int *data_int;
				unsigned long *data_long;
				double range;
				double offset;

				if ( sizeof(int) == 4 )
				{
					data_int = new unsigned int[numPixels2];
				}
				else if ( sizeof(long) == 4 )
				{
					data_long = new unsigned long[numPixels2];
				}
				else
				{
					cout << "ERROR. In edf::write_data() - Neither 'unsigned int' nor 'unsigned long int' are of size 4 byte." << endl;
					return 3;
				}

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( std::pow( 256.0, sizeDataType1 ) - 1.0 );
					offset = get_ScalingMin();
				}
				else
				{
					range  = 1;
					offset = 0;
				}
				for (unsigned long i = 0; i < numPixels2; i++)
				{
					if ( range != 0 )
					{
						if ( sizeof(int) == 4 )
						{
							data_int[i] = (unsigned int )floor( (mat[i]-offset)/range+epsilon);
						}
						else
						{
							data_long[i] = (unsigned long )floor( (mat[i]-offset)/range+epsilon);
						}
					}
					else
					{
						if ( sizeof(int) == 4 )
						{
							data_int[i] = 0;
						}
						else
						{
							data_long[i] = 0;
						}
					}
				}		

				if ( reverse_ByteOrder )
				{
					if ( sizeof(int) == 4 )
					{
						revByteOrder( data_int, numPixels2 );
					}
					else
					{
						revByteOrder( data_long, numPixels2 );
					}
				}

				if ( sizeof(int) == 4 )
				{
					fout.write( (char *)data_int, numBytes );
				}
				else
				{
					fout.write( (char *)data_long, numBytes );
				}
				if ( !fout )
				{
					cout << "WARNING. Problems appeared writing to " << get_FilenameOut() << endl;
				}

				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numPixels1 << complex_str << "unsigned longs written to "
						<< get_FilenameOut() << " (" << un << "scaled)" << endl;
				}

				if ( sizeof(int) == 4 )
				{
					delete []data_int;
				}
				else
				{
					delete []data_long;
				}
			}
			break;

			case DT_NOT_SPECIFIED:
			{

			}
			break;

			case DT_FAMELESS:
			{

			}
			break;

			default:
			{
				// The user should never see this message.
				cerr << "ERROR. In edf::get_datatype_str() - ";
				cerr <<	"Invalid data type '" << get_DataType() << "'." << endl;
				set_DataType( DT_NOT_SPECIFIED );
				retval = 1;
			}
			break;
		}

		return retval;
	}

	//--------------------------------------------------------------------- read
	//
	//
	int edf::read( double *&mat, bool flipByteOrder )
	{
		return read( mat, get_FilenameIn(), flipByteOrder );
	}

	//--------------------------------------------------------------------- read
	//
	//
	int edf::read( double *&mat, string filename, bool flipByteOrder )
	{
		int retval;
		
		retval = read_header( filename );
		if ( retval )
		{
			return 1;
		}
		
		//cout << "1==> dimensions: " << get_NumberOfDimensions() << endl;
		retval = read_data( mat, filename, flipByteOrder );
		//cout << "2==> dimensions: " << get_NumberOfDimensions() << endl;
		if ( retval )
		{
			return 2;
		}

		return 0;
	}

	//-------------------------------------------------------------------- write
	//
	//
	int edf::write( double *mat, bool flipByteOrder )
	{
		return write( mat, get_FilenameOut(), flipByteOrder );
	}

	//-------------------------------------------------------------------- write
	//
	//
	int edf::write( double *mat, string filename, bool flipByteOrder )
	{
		int retval;

		retval = write_header( filename );
		if ( retval )
		{
			return 1;
		}
	
		retval = write_data( mat, filename, flipByteOrder );
		if ( retval )
		{
			return 2;
		}

		return 0;
	}


	//--------------------------------------------------------------- parse_line
	//
	//
	int edf::parse_line( string input, C_Tagpair &tagpair ) const
	{
		string::size_type pos0;
		string::size_type pos1;
		string::size_type pos2;
		string::size_type pos3;
		string::size_type pos4;
		string::size_type pos5;
		string::size_type length_end;
		string::size_type length_assign;
		string::size_type length_separator;
		C_Keyword keyword;
		string keyname;
		string value;
		string comment;

		length_end = string(HEADER_TOKEN_END).length();
		length_separator = string(HEADER_TOKEN_SEPARATOR).length();
		length_assign = string(HEADER_TOKEN_ASSIGN).length();

		tagpair.reset();

		if ( input.length() == 0 )
		{
			// The line is empty. Skip this line.
			return 0;
		}

		if ( input.find( HEADER_TOKEN_END ) == input.length()-length_end )
		{
			// The lines indicates the end of the header.
			return 1;	
		}

		pos0 = input.rfind( HEADER_TOKEN_SEPARATOR );
		if ( pos0 == string::npos )
		{
			// Missing default separator --> try alternative one
			pos0 = input.rfind( HEADER_TOKEN_SEPARATOR_ALT );
			if ( pos0 == string::npos )
			{
				// Alternative separator also missing --> syntax error.
				return 2;
			}
			length_separator = string(HEADER_TOKEN_SEPARATOR_ALT).length();
		}
		
		comment = input.substr( pos0 + length_separator );
		
		// The comparison '>=' must not be replaced by '==', since the separator
		// may start with whitespaces. (In the current version it is " ;".)
		if ( input.find_first_not_of( " \t" ) >= pos0 )
		{
			// There are only whitespace characters before the separator.
			// Ignore this line.
			return 0;
		}

		pos1 = input.find( HEADER_TOKEN_ASSIGN );
		if ( pos1 == string::npos )
		{
			// Missing assign token --> syntax error.
			return 3;	
		}
		
		keyname = input.substr( 0, pos1 );
		pos2 = keyname.find_first_not_of( " \t" );
		if ( pos2 == string::npos )
		{
			// Missing keyword before assign token --> syntax error.
			return 4;	
		}

		pos3 = keyname.find_first_of( " \t" );
		keyname = keyname.substr( pos2, pos3-pos2 );

		value = input.substr( pos1+length_assign, pos0-pos1-length_assign );
		pos4 = value.find_first_not_of( " \t" );
		pos5 = value.find_last_not_of( " \t" );
		
		if ( pos4 != string::npos )
		{
			// There is not only whitespace between assign token and separator.
			// This means that the given keyword has got a value. The value is
			// given by the part between the first and last non whitespace.
			// (All whitespace between is considered being part of the value.
			
			value = value.substr( pos4, pos5-pos4+1 );
		}
		
		keyword[0] = keyname;
		for ( C_AllKeywords::const_iterator i = p_allKeywords.begin();
			  i != p_allKeywords.end(); i++ )
		{
			if ( keyname == *i )
			{
				keyword = *i;
				break;
			}
		}
		tagpair.set( keyword, value, comment );

		return 0;
	}

	//---------------------------------------------------------- interpret_tags
	//
	//
	int edf::interpret_tags()
	{
		tags_t tags;
		bool newflag;
		bool matlibflag; // Indicates that edf file had been written using matlib.

		//set_DataType( DT_NOT_SPECIFIED );
		//set_ByteOrderIn( BO_NOT_SPECIFIED );
		//set_RealFormat( REAL_FORMAT_NOT_SPECIFIED );
		set_ScaleIn( false );	// Will be set to 'true' if RWTHscaling is defined in edf header.
								// RWTHscaling is true for tagname 'TAG_TUD_SCALING_MAX' is defined.
								// Question by JP: would it be wise to have a tag 'RWTHscaling'?

		newflag = false;
		matlibflag = false;
		for ( int i=0; i<p_tags_loaded.size(); i++ )
		{
			C_Tagpair tag;
			string mainKey;
			string value_str;
			bool found;
			int direction=0;
			
			tag = p_tags_loaded.get( i );
			mainKey = tag.get_mainKey();

			for ( C_Keyword::const_iterator i = TAG_DIM.begin();
				  i != TAG_DIM.end(); i++ )
			{
				if ( mainKey.find( *i ) == 0 )
				{
					// Direction in {0, 1, 2, ...}
					// Keyword in {Dim_1, Dim_2, Dim_3, ...}.
					// So do not forget the minus 1.
					direction = toInt( mainKey.substr( i->length() ) ) - 1;
					mainKey = *i;
					//cout << "DEBUGINFO. In edf::interpret_tags() - ";
					//cout << "mainKey = '" << mainKey << "'" << endl;
					//cout << "\t" << "direction = " << direction << endl;
				}
			}

			found = true;
			if		( mainKey == TAG_HEADER_ID )
			{
				set_HeaderID( tag.to_string() );
			}
			else if ( mainKey == TAG_TUD_HEADER_SIZE )
			{
				// identifies it as standard idendifier but does notthing! CS
			}
			else if ( mainKey == TAG_TUD_HEADER_LINES )
			{
				// identifies it as standard idendifier but does notthing! CS
			}
			else if ( mainKey == TAG_SIZE )
			{
				set_Size( tag.to_ulint() );
			}
			else if ( mainKey == TAG_DIM )
			{
				//cout << "DEBUGINFO. In edf::interpret_tags() - ";
				//cout << "A number of dimensions = " << get_NumberOfDimensions()  << endl;
				set_Dimension( direction, tag.to_ulint() );
				//cout << "DEBUGINFO. In edf::interpret_tags() - ";
				//cout << "B number of dimensions = " << get_NumberOfDimensions()  << endl;
			}
			else if ( mainKey == TAG_BYTE_ORDER )
			{
				set_ByteOrderIn( tag.to_string() );
			}
			else if ( mainKey == TAG_DATA_TYPE )
			{
				set_DataType( tag.to_string() );
			}
			else if ( mainKey == TAG_REAL_FORMAT )
			{
				set_RealFormat( tag.to_string() );
			}
			else if ( mainKey == TAG_TUD_SCALING_MIN )
			{
				set_ScalingMin( tag.to_double() );
			}
			else if ( mainKey == TAG_TUD_SCALING_MAX )  							
			{ 
				set_ScalingMax( tag.to_double() );
				set_ScaleIn( true );
			}
			else if ( mainKey == TAG_TUD_FILE_TYPE )
			{
				set_FileType( tag.to_string() );
				matlibflag = true;
			}
			else if ( mainKey == TAG_TUD_DATE_ORIG )
			{
				set_DateOrig( tag.to_string() );
				newflag = true;
			}
			else if ( mainKey == TAG_TUD_DATE_PREV )
			{
				// Identifies it as built-in keyword but does nothing. JP
				newflag = true;
			}
			else if ( mainKey == TAG_TUD_DATE_CURR )
			{
				set_DatePrev( tag.to_string() );
				newflag = true;
			}
			else if ( mainKey == TAG_TUD_FILENAME_ORIG	)
			{
				set_FilenameOrig( tag.to_string() );
			}
			else if ( mainKey == TAG_TUD_FILENAME_PREV	)
			{
				// Identifies it as built-in keyword but does nothing. JP
			}
			else if ( mainKey == TAG_TUD_FILENAME_CURR	)
			{
				// Identifies it as built-in keyword but does nothing. JP
			}
			else if ( mainKey == TAG_IMAGE )
			{
				set_ImageNumber( tag.to_int() );
			}
			else if ( mainKey == TAG_VERSION_NUMBER )
			{
				set_VersionNumber( tag.to_string() );
			}
			else
			{
				found = false;
			}
			
			if ( found )
			{
				p_tags_builtin.set( tag );
			}
			else
			{
				p_tags_userdef.set( tag );
			}
		}

		// The newflag is set in oder to distinguish between the old (errornous)
		// and the new (correct) byteorder specification. If old then byteorder
		// must be flipped in case of double or float.
		set_NewFlag( newflag && matlibflag );
		set_MatlibFlag( matlibflag );

		return (int )!check_consistency();
	}

	//----------------------------------------------------- calculate_headersize
	//
	//
	unsigned long edf::calculate_headersize( text_t text, int numCharsEndl, int &numSpacePadding ) const
	{
	/*
		Some comments on the algorithm used to determine the header size including meta information.

		Definition:
		header: ASCII encoded part at beginning of edf image starting with '{\n' ending with '}\n'.
		headersize: size of the header measured in bytes
		meta part: all characters of the header describing the headersize
		space part: all characters of the header needed for filling to next multiple of HEADER_BYTES_MODULO
		normal part: all characters of the header that do not belong to the meta part nor to the space part

		Symbols:
		N: headersize
		a: size of the normal part
		b: size of the meta part
		c: size of the space part

		Functions:
		log: logarithm with respect to basis 10
		len: len(x) = floor( log(n) ) + 1, i.e. number of digits for x in decimal representation
		pad: pad(x) = ((x-1)/2^10+1)*2^10, i.e. smallest number n*2^10 with x<=n*2^10 ['/' means integer division]

		Hence:
		N = a + b + c
		N = pad( a + b )
		len(N) = b
		len( pad(a+b) ) = b
		floor( log( pad(a+b) ) ) + 1 = b
		floor( log( ((a+b-1)/2^10+1)*2^10 ) ) + 1 = b   (*)

		Given the value of a, we need to find a value of b solving (*).
		In order to solve (*) numerically, I define following sequence:
		b(0) = 0
		b(n+1) = floor( log( ((a+b(n)-1)/HEADER_BYTES_MODULO+1)*HEADER_BYTES_MODULO ) ) + 1   (**)
		The limit of this sequence b(n) is a solution of (*).

		There are 3 important facts about (**):
		1. b(n) converges, i.e. the limit b exists.
		2. b is element of {b(n), n>=0}, i.e. the limit will be reached by the series.
		3. b is element of {b(n), n in {1, 2}}, i.e. the limit will be reached with only 2 or 3 iterations.
		4. The case of 3 iterations only appears for a >= (10^10-2^10) - (10-2) = (5^10-1)*2^10 - (10-2) = 9999998968
		   which is too big to be stored in a <long int> variable...

		Therefore using (**) will retrieve the required headersize with exactly 2 iterations.

		by JP in May 2007
	*/

		long int a;
		long int b;
		long int c;
		long int b_old;
		long int headersize;
		int count;

		a = 0;
		for ( text_t::size_type i=0; i<text.size(); i++ )
		{
			a += (int )text[i].length() + numCharsEndl;
		}
		a -= HEADER_WIDTH_VALUE; // new
		//cout << "DEBUGINFO: In edf::calculate_headersize() - before: a=" << a << endl;
		count = 0;
		b = 0;
		do
		{
			//long int N;

			count++;
			b_old = b;
	;
			const long int N = ( (a+b-1) / HEADER_BYTES_MODULO + 1 ) * HEADER_BYTES_MODULO;
			b = (long int )floor( log10( (long double) N ) ) + 1;
			c = N - a - b;

			headersize = N;
		} while ( b != b_old );

		//cout << "DEBUGINFO: In edf::calculate_headersize() - after: a=" << a << endl;
		//cout << "DEBUGINFO: In edf::calculate_headersize() - after: b=" << b << endl;
		//cout << "DEBUGINFO: In edf::calculate_headersize() - after: c=" << c << endl;

		if ( b < HEADER_WIDTH_VALUE )
		{
			numSpacePadding = (int)(c - HEADER_WIDTH_VALUE + b); // new
			
			// 01.12.2007 3:00, bug fix at bw4
			while ( numSpacePadding < 0 )
			{
				numSpacePadding += HEADER_BYTES_MODULO;
			}
		}
		else
		{
			numSpacePadding = (int)c;
		}
		
		return headersize;
	}

	//------------------------------------------------------- calculate_datasize
	//
	//
	unsigned long edf::calculate_datasize() const
	{
		unsigned long datasize;
		
		if ( get_NumberOfDimensions() == 0 )
		{
			datasize = 0;
		}
		else
		{
			datasize = get_SizeOfDataType();
			for ( int i=0; i<get_NumberOfDimensions(); i++ )
			{
				datasize *= get_Dimension(i);
			}
		}
		
		return datasize;
	}


	//--------------------------------------------------------------------- checkConsistency
	//
	//
	bool edf::check_consistency( bool readonly )
	{
		bool consistencyflag = true;

		if ( get_ByteOrderIn() == BO_NOT_SPECIFIED )
		{
			cerr << "WARNING. No ByteOrder specified. Assuming '" << HIGH_BYTE_FIRST_STR << "'." << endl;
			if ( !readonly )
			{
				set_ByteOrderIn( BO_HIGH_BYTE_FIRST );
			}
		}

		if ( get_ByteOrderIn() == BO_FAMELESS )
		{
			cerr << "WARNING. Unknown ByteOrder. Assuming '" << HIGH_BYTE_FIRST_STR << "'." << endl;
			if ( !readonly )
			{
				set_ByteOrderIn( BO_HIGH_BYTE_FIRST );
			}
		}

		if ( get_DataType() == DT_NOT_SPECIFIED )
		{
			cerr << "WARNING. DataType not specified. Assuming '" << UNSIGNED_SHORT_STR << "'." << endl;
			if ( !readonly )
			{
				set_DataType( DT_UNSIGNED_SHORT );
			}
		}

		if ( get_DataType() == DT_FAMELESS )
		{
			cerr << "WARNING. Unknown datatype. Assuming '" << UNSIGNED_SHORT_STR << "'." << endl;
			if ( !readonly )
			{
				set_DataType( DT_UNSIGNED_SHORT );
			}
		}

		if ( ( get_DataType() < datatype_min ) || ( get_DataType() > datatype_max ) )
		{
			// The user should never see this message.
			cerr << "ERROR! Bug in edf::getSizeOfDataType(). ";
			cerr <<	"Invalid data type '" << get_DataType() << "'." << endl;
			if ( !readonly )
			{
				set_DataType( DT_NOT_SPECIFIED );
			}
		}
		
		if ( is_datatype_real() )
		{
			if ( get_RealFormat() == REAL_FORMAT_NOT_SPECIFIED )
			{
				cerr << "WARNING. RealFormat not specified. Assuming '" << REAL_FORMAT_IEEE_STR << "'." << endl;
				if ( !readonly )
				{
					set_RealFormat( REAL_FORMAT_IEEE );
				}
			}

			if ( get_RealFormat() == REAL_FORMAT_FAMELESS )
			{
				cerr << "WARNING. Unknown RealFormat. Assuming '" << REAL_FORMAT_IEEE_STR << "'." << endl;
				if ( !readonly )
				{
					set_RealFormat( REAL_FORMAT_IEEE );
				}
			}
		}
		
		// Eliminated beacuse a zero data block is allowed. (JP)
		/*
		if (!get_Size())
		{
			cerr << "WARNING. Data size not specified. Assuming Dim1*Dim2*sizeof(datatype)!" << endl;
			set_Size(getDim_x() * get_Dim_y() * get_SizeOfDataType());
		}
		*/
		
		if ( get_Size() != calculate_datasize() )
		{
			cerr << "ERROR. Size not compatible with image dimensions!" << endl;
			cerr << "Size of datatype: " << get_datatype_string(get_DataType())
				 << " (" << get_SizeOfDataType() << " bytes)" << endl;
			cerr << "Size given by header attribute:  " << get_Size() << endl;
			cerr << "Size calculated from dimensions: " << calculate_datasize() << endl;
			consistencyflag = false;
		}

		// Eliminated because TomoAngleRangeRWTH is not a built-in tag any more.
		/*
		if (getTomoAngleRangeRWTH() > 0)
		{
			if (fabs(getTomoAngleRangeRWTH()*180. - (get_Dim_y()-1)*getYStepRWTH()) > 1.e-3)
				cerr << "WARNING. Angular scaling not correct!" << endl;
		}
		*/

		return consistencyflag;
	}

	//-------------------------------------------------------------- writeTagStr
	//
	//
	int edf::writeTagStr( text_t &text, C_Keyword keyword, string value, string comment )
	{
		string mainKey;
		
		mainKey = keyword.get_main();
		text.push_back( mainKey + spacepadding( mainKey, HEADER_WIDTH_KEYWORD ) + HEADER_TOKEN_ASSIGN +
						value + spacepadding( value, HEADER_WIDTH_VALUE ) + HEADER_TOKEN_SEPARATOR +
						comment );

		return 0;
	}

	//--------------------------------------------------------- writeTagStrIndex
	//
	//
	int edf::writeTagStrIndex( text_t &text, text_t::size_type index, C_Keyword keyword, string value, string comment )
	{
		string mainKey;

		if ( index >= text.size() )
		{
			return 1;
		}

		mainKey = keyword.get_main();
		text[index] = ( mainKey + spacepadding( mainKey, HEADER_WIDTH_KEYWORD ) + HEADER_TOKEN_ASSIGN +
						value + spacepadding( value, HEADER_WIDTH_VALUE ) + HEADER_TOKEN_SEPARATOR +
						comment );

		return 0;
	}

	//-------------------------------------------------------------------- print
	//
	//
	string edf::print() const
	{
		ostringstream osst;

		osst << "Built-in tags:" << endl;
		osst << p_tags_builtin << endl;
		osst << "User-defined tags:" << endl;
		osst << p_tags_userdef << endl;
		osst << "Tags as they have been loaded:" << endl;
		osst << p_tags_loaded << endl;
		
		return osst.str();
	}
	
	//==========================================================================
	//========================================================= Global functions

	//-------------------------------------------------------------get_datatypes
	//
	//
	string get_datatypes()
	{
		ostringstream osst;

		osst << " Float= "				<< DT_FLOAT;
		osst << " Double="				<< DT_DOUBLE;

		osst << " SignedByte="			<< DT_BYTE;
		osst << " SignedShort=" 		<< DT_SHORT;
		osst << " SignedInteger="		<< DT_INTEGER;
		osst << " SignedLong="			<< DT_LONG;
															
		osst << " UnsignedByte="		<< DT_UNSIGNED_BYTE;
		osst << " UnsignedShort="		<< DT_UNSIGNED_SHORT;
		osst << " UnsignedInteger=" 	<< DT_UNSIGNED_INTEGER;
		osst << " UnsignedLong="		<< DT_UNSIGNED_LONG;
															
		osst << " ComplexFloat="		<< DT_COMPLEX_FLOAT;
		osst << " ComplexDouble="		<< DT_COMPLEX_DOUBLE;
		osst << " ComplexByte=" 		<< DT_COMPLEX_BYTE;
		osst << " ComplexShortInteger=" << DT_COMPLEX_SHORT;
		osst << " ComplexInteger="		<< DT_COMPLEX_INTEGER;
		osst << " ComplexLongInteger="	<< DT_COMPLEX_LONG;

		return osst.str();
	}
	
	//------------------------------------------------------ get_datatype_string
	//
	//
	string get_datatype_string( ns_edf::datatype_t type )
	{
		string retstring;
		
		switch ( type )
		{
			case DT_FLOAT:				retstring = "DT_FLOAT"; break;
			case DT_DOUBLE: 			retstring = "DT_DOUBLE"; break;

			case DT_BYTE:				retstring = "DT_BYTE"; break;
			case DT_SHORT:				retstring = "DT_SHORT"; break;
			case DT_INTEGER:			retstring = "DT_INTEGER"; break;
			case DT_LONG:				retstring = "DT_LONG"; break;

			case DT_UNSIGNED_BYTE:		retstring = "DT_UNSIGNED_BYTE"; break;
			case DT_UNSIGNED_SHORT:		retstring = "DT_UNSIGNED_SHORT"; break;
			case DT_UNSIGNED_INTEGER:	retstring = "DT_UNSIGNED_INTEGER"; break;
			case DT_UNSIGNED_LONG:		retstring = "DT_UNSIGNED_LONG"; break;

			case DT_COMPLEX_FLOAT:		retstring = "DT_COMPLEX_FLOAT"; break;
			case DT_COMPLEX_DOUBLE:		retstring = "DT_COMPLEX_DOUBLE"; break;
			case DT_COMPLEX_BYTE:		retstring = "DT_COMPLEX_BYTE"; break;
			case DT_COMPLEX_SHORT:		retstring = "DT_COMPLEX_SHORT"; break;
			case DT_COMPLEX_INTEGER:	retstring = "DT_COMPLEX_INTEGER"; break;
			case DT_COMPLEX_LONG:		retstring = "DT_COMPLEX_LONG"; break;

			case DT_NOT_SPECIFIED:		retstring = "DT_NOT_SPECIFIED"; break;
			case DT_FAMELESS:			retstring = "DT_FAMELESS"; break;

			default: retstring = "unknown";
		}

		return retstring;
	}

	//---------------------------------------------------------------- operator<<
	//
	//
	ostream &operator<<( ostream &os, const edf &file )
	{
		return os << file.print();
	}

}

//==============================================================================
//==================================================================== OLD STUFF
//==============================================================================
/*
	int edf::read(double *mat, bool flipByteOrder)
	{
		string filename("");
 		ifstream fin;
		bool wait;
		unsigned long length;
		unsigned long index;
		int retval;

		retval = 0;
		wait = true;
		filename = getFilename();
		fin.open( filename.c_str() );
		while ( (!fin || fin.bad() || !fin.good()) && wait )
		{
			ns_commandline::commandline cline;

			cline << "ERROR. Could not open file '" + filename + "'. Wait?";
			cline.set_hint( "0/1" );
			cline.set_default( 1 );
			cline >> wait;

			fin.open( filename.c_str() );
			matlibutil::usleep( 200 * 1000 );
		}

		if ( !fin || fin.bad() || !fin.good() )
		{
			cerr << "ERROR. Could not open file '" + filename + "'." << endl;
			return 1;
		}
		else
		{
			readHeader(filename);
			char dummy;

			while ( !fin.eof() && ( fin.get() != '}' ) );
			while ( !fin.eof() && ( fin.get() != 10  ) );

			switch ( getDataType() )
			{
				case DT_DOUBLE:
				{
					if ( (ny != get_Dim_y()) || (nx != get_Dim_x()) )
					{
						ny = get_Dim_y();
						nx = get_Dim_x();
						delete mat;
						mat = new double[get_Dim_y()*get_Dim_x()];
					}

					length = (unsigned long) fread(mat, sizeof(double), get_Dim_x()*get_Dim_y(), f);

					if ( getVerbosity() )
						cout << length << " doubles read" << endl;

					if (length != get_Dim_x()*get_Dim_y())
						cerr << "WARNING. There is trouble with the file size!" << endl;

					if (BigEndianByteOrder())
						revByteOrder(mat, get_Dim_y()*get_Dim_x());
					if (flipByteOrder)
						revByteOrder(mat, get_Dim_y()*get_Dim_x());
				}
				break;

				case DT_FLOAT:
				{
					float *data;

					data = new float[get_Dim_x()*get_Dim_y()];
					length = (unsigned long) fread(data, sizeof(float), get_Dim_x()*get_Dim_y(), f);

					if ( getVerbosity() )
						cout << length << " floats read" << endl;

					if (length != get_Dim_x()*get_Dim_y())
						cerr << "WARNING. There is trouble with the file size!" << endl;

					if ((ny != get_Dim_y()) || (nx != get_Dim_x()))
					{
						ny = get_Dim_y();
						nx = get_Dim_x();
						delete mat;
						mat = new double[get_Dim_y()*get_Dim_x()];
					}

					for (unsigned int i = 0; i < ny; i++)
					{
						for (unsigned int j = 0; j < nx; j++)
						{
							index = (unsigned long )i * (unsigned long )nx + (unsigned long )j;

							if (BigEndianByteOrder())
								revByteOrder(&(data[index]));
							if (flipByteOrder)
								revByteOrder(&(data[index]));
							mat[index] = data[index];
						}
					}

					delete []data;		
				}
				break;

				case DT_SHORT:
				{
					short *data;
					double range;
					double offset;
					double shift;

					data = new  short[get_Dim_x()*get_Dim_y()];
					length = (unsigned long) fread(data, sizeof(short), get_Dim_x()*get_Dim_y(), f);

					if ( getVerbosity() )
						cout << length << " shorts read " << endl;

					if (length != get_Dim_x()*get_Dim_y())
						cerr << "WARNING. There is trouble with the file size!" << endl;
					if (flipByteOrder)
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ((ny != get_Dim_y()) || (nx != get_Dim_x()))
					{
						ny = get_Dim_y();
						nx = get_Dim_x();
						delete mat;
						mat = new double[get_Dim_y()*get_Dim_x()];
					}

					if ( get_ScaleIn() )
					{
						range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 - 1.0 );
						offset = get_ScalingMin();
						shift  = 65536.0 / 2.0;
					}
					else
					{
						range  = 1;
						offset = 0;
						shift  = 0;
					}

					for (unsigned int i = 0; i < ny; i++)
					{
						for (unsigned int j = 0; j < nx; j++)
						{
							index = (unsigned long )i * (unsigned long )nx + (unsigned long )j;
							mat[index] = ( (double )data[index] + shift ) * range + offset;
						}
					}
					delete []data;
				}
				break;

				case DT_UNSIGNED_SHORT:
				{
					unsigned short *data;
					double range;
					double offset;

					data = new unsigned short[get_Dim_x()*get_Dim_y()];
					length = (unsigned long) fread(data, sizeof(unsigned short), get_Dim_x()*get_Dim_y(), f);

					if ( getVerbosity() )
						cout << length << " (unsigned) shorts read " << endl;

					if (length != get_Dim_x()*get_Dim_y())
						cerr << "WARNING. There is trouble with the file size!" << endl;
					if (flipByteOrder)
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ((ny != get_Dim_y()) || (nx != get_Dim_x()))
					{
						ny = get_Dim_y();
						nx = get_Dim_x();
						delete mat;
						mat = new double[get_Dim_y()*get_Dim_x()];
					}

					if ( get_ScaleIn() )
					{
						range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 - 1.0 );
						offset = get_ScalingMin();
					}
					else
					{
						range  = 1;
						offset = 0;
					}

					for (unsigned int i = 0; i < ny; i++)
					{
						for (unsigned int j = 0; j < nx; j++)
						{
							index = (unsigned long )i * (unsigned long )nx + (unsigned long )j;
							mat[index] = (double )data[index] * range + offset;
						}
					}
					delete []data;
				}
				break;

				case DT_LONG:
				{
					long *data;
					double range;
					double offset;
					double shift;

					data = new int long[get_Dim_x()*get_Dim_y()];
					length = (unsigned long) fread(data, sizeof(long), get_Dim_x()*get_Dim_y(), f);

					if ( getVerbosity() )
						cout << length << " longs read " << endl;

					if (length != get_Dim_x()*get_Dim_y())
						cerr << "WARNING. There is trouble with the file size!" << endl;
					if (flipByteOrder)
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ((ny != get_Dim_y()) || (nx != get_Dim_x()))
					{
						ny = get_Dim_y();
						nx = get_Dim_x();
						delete mat;
						mat = new double[get_Dim_y()*get_Dim_x()];
					}

					if ( get_ScaleIn() )
					{
						range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 * 65536.0 - 1.0 );
						offset = get_ScalingMin();
						shift  = ( 65536.0 * 65536.0 ) / 2.0;
					}
					else
					{
						range  = 1;
						offset = 0;
						shift  = 0;
					}

					for (unsigned int i = 0; i < ny; i++)
					{
						for (unsigned int j = 0; j < nx; j++)
						{
							index = (unsigned long )i * (unsigned long )nx + (unsigned long )j;
							mat[index] = ( (double )data[index] + shift ) * range + offset;
						}
					}
					delete []data;
				}
				break;

				case DT_UNSIGNED_LONG:
				{
					unsigned long *data;
					double range;
					double offset;

					data = new unsigned long[get_Dim_x()*get_Dim_y()];
					length = (unsigned long) fread(data, sizeof(unsigned long), get_Dim_x()*get_Dim_y(), f);

					if ( getVerbosity() )
						cout << length << " (unsigned) longs read " << endl;

					if (length != get_Dim_x()*get_Dim_y())
						cerr << "WARNING. There is trouble with the file size!" << endl;
					if (flipByteOrder)
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
						revByteOrder(data, length);
					if ((ny != get_Dim_y()) || (nx != get_Dim_x()))
					{
						ny = get_Dim_y();
						nx = get_Dim_x();
						delete mat;
						mat = new double[get_Dim_y()*get_Dim_x()];
					}

					if ( get_ScaleIn() )
					{
						range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 * 65536.0 - 1.0 );
						offset = get_ScalingMin();
					}
					else
					{
						range  = 1;
						offset = 0;
					}

					for (unsigned int i = 0; i < ny; i++)
					{
						for (unsigned int j = 0; j < nx; j++)
						{
							index = (unsigned long )i * (unsigned long )nx + (unsigned long )j;
							mat[index] = (double )data[index] * range + offset;
						}
					}
					delete []data;
				}
				break;

				case DT_NOT_SPECIFIED:
				{

				}
				break;

				case DT_FAMELESS:
				{

				}
				break;

				default:
				{
					// The user should never see this message.
					cerr << "ERROR! Bug in edf::get_datatype_str(). ";
					cerr <<	"Invalid data type '" << getDataType() << "'." << endl;
					setDataType( DT_NOT_SPECIFIED );
					retval = 1;
				}
				break;
			}
		}
		fclose(f);
		return retval;
	}

	int edf::writeTo(unsigned long int ny, unsigned long int nx, double *mat, string fn)
	{
		setFilename( fn );
		return(write(ny, nx, mat));
	}

	void edf::setParams(unsigned long int ny, unsigned long int nx)
	{
		setDim_1(nx);
		setDim_2(ny);
		setSize(ny*nx*getSizeOfDataType());
	}

	int edf::write(unsigned long int ny, unsigned long int nx, double *mat, bool flipByteOrder)
	{
		string filename("");
 		ofstream fout;
		text_t t;
		int retval;
		int headersize;
		int numCharsEndl;
		int numSpacePadding;
		fpos_t filepos;
		unsigned long numBytesWritten;
		unsigned long index;
		double epsilon; 	// Added to double values before they are converted to integers
							// in order to ensure the scaled maximum is 2^n-1 instead of 2^n-2

		epsilon = 0.1;
		retval = 0;

		filename = getFilename();
		fout.open( filename.c_str() );
		if ( fout || fout.bad() || !fout.good() )
		{
			cerr << "ERROR. Could not open file '" << filename << "'. Data not written." << endl;
			cin.ignore();
			return 1;
		}

		setParams(ny, nx);
		
		// Retrieve the platform depentent number of bytes used to represent a newline.		
		fprintf(f, "{\n");
		fgetpos(f, &filepos);
		#ifdef linux
			numCharsEndl = filepos.__pos - 1; // -1 because of the character '{'
		#else
			numCharsEndl = (int)filepos - 1; // -1 because of the character '{'
		#endif

		// Write the header in datastructure text_t in order to have random access.
		// Random access is needed to fill in the headersize value after it will have
		// been calculated. But for calculating headersize, we need to write the
		// header first... (See function calculateHeadersize().)
		t.clear();
		t.push_back( "{" );
		writeTag( t, TAG_HEADER_ID, 			getHeaderID(), 				isHeaderID_defined()	 			);
		writeTag( t, TAG_HEADER_SIZE_TUD,		"",							isHeaderSizeTUD_defined()			);
		writeTag( t, TAG_HEADER_LINES_TUD,		"",							isHeaderLinesTUD_defined()			);
		writeTag( t, TAG_IMAGE, 				getImage(),					isImage_defined()   				);
		writeTag( t, TAG_BYTE_ORDER, 			getByteOrder_str(),			isByteOrder_defined()  				);
		writeTag( t, TAG_DATA_TYPE, 			getDataType_str(),			isDataType_defined()   				);
		writeTag( t, TAG_DIM_1, 				get_Dim_x(),					isDim_1_defined()    				);
		writeTag( t, TAG_DIM_2, 				get_Dim_y(),					isDim_2_defined()		   			);
		writeTag( t, TAG_SIZE,  				getSize(),					isSize_defined()		   			);
		writeTag( t, TAG_TITLE, 				getTitle(),					isTitle_defined()					);
		writeTag( t, TAG_PREFIX, 				getPrefix(), 				isPrefix_defined()					);
		writeTag( t, TAG_SUFFIX, 				getSuffix(), 				isSuffix_defined()					);
		writeTag( t, TAG_VERSION_RWTH,			getVersionRWTH(),			isVersionRWTH_defined()				);
		writeTag( t, TAG_FILE_TYPE_RWTH,		getFileTypeRWTH(),			isFileTypeRWTH_defined()			);
		writeTag( t, TAG_DATE_TUD,				getDateTUD(),				isDateTUD_defined()					);
		writeTag( t, TAG_COUNT_TIME,			getCount_time(),			isCount_time_defined()	    		);
		writeTag( t, TAG_ACCTR_DATE_TUD,		getAcctrDateTUD(),			isAcctrDateTUD_defined()			);
		writeTag( t, TAG_ACCTR_ENVIRON_TUD,		getAcctrEnvironTUD(),		isAcctrEnvironTUD_defined()			);
		writeTag( t, TAG_ACCTR_SESSION_TUD, 	getAcctrSessionTUD(),		isAcctrSessionTUD_defined()			);
		writeTag( t, TAG_ACCTR_LOGFILE_TUD, 	getAcctrLogfileTUD(),		isAcctrLogfileTUD_defined()			);
		writeTag( t, TAG_DETECTOR_NAME_TUD,		getDetectorNameTUD(),		isDetectorNameTUD_defined()			);
		writeTag( t, TAG_X_BINNING_TUD, 		getXBinningTUD(),			isXBinningTUD_defined()				);
		writeTag( t, TAG_Y_BINNING_TUD, 		getYBinningTUD(),			isYBinningTUD_defined()				);
		writeTag( t, TAG_TEMPERATURE_TUD,		getTemperatureTUD(),		isTemperatureTUD_defined()			);
		writeTag( t, TAG_SCAN_NO, 				getScan_no(),				isScan_no_defined()	    			);
		writeTag( t, TAG_POINT_NO,  			getPoint_no(),				isPoint_no_defined()	    		);
		writeTag( t, TAG_PRESET, 				getPreset(),				isPreset_defined() 	    			);
		writeTag( t, TAG_COUNTER_MNE,			getCounter_mne(),			isCounter_mne_defined()				);
		writeTag( t, TAG_COUNTER_POS, 			getCounter_pos(),			isCounter_pos_defined()				);
		writeTag( t, TAG_MOTOR_MNE, 			getMotor_mne(),				isMotor_mne_defined()				);
		writeTag( t, TAG_MOTOR_POS, 			getMotor_pos(),				isMotor_pos_defined()				);
		writeTag( t, TAG_SR_CURRENT,			getSRCur(), 				isSRCur_defined()					);
		writeTag( t, TAG_RUN,					getRun(),					isRun_defined()						);
		writeTag( t, TAG_MAX_VAL_RWTH,			get_ScalingMax(),			isMaxValRWTH_defined() 				);
		writeTag( t, TAG_MIN_VAL_RWTH,			get_ScalingMin(),			isMinValRWTH_defined() 				);
		writeTag( t, TAG_X_STEP_RWTH,			getXStepRWTH(), 			isXStepRWTH_defined()				);
		writeTag( t, TAG_Y_STEP_RWTH,			getYStepRWTH(), 			isYStepRWTH_defined() 				);
		writeTag( t, TAG_ZERO_PAD_RWTH, 		getZeroPadRWTH(),			isZeroPadRWTH_defined()				);
		writeTag( t, TAG_ROT_AXIS_RWTH, 		getRotAxisRWTH(),			isRotAxisRWTH_defined()				);
		writeTag( t, TAG_ROT_AXIS_X_POS_RWTH,	getRotAxisXPosRWTH(),		isRotAxisXPosRWTH_defined()			);
		writeTag( t, TAG_ROT_AXIS_X_POS_RWTH,	getRotAxisYPosRWTH(),		isRotAxisYPosRWTH_defined()			);
		writeTag( t, TAG_ROT_ANGLE_RWTH,		getRotAngleRWTH(),			isRotAngleRWTH_defined()			);
		writeTag( t, TAG_TOMO_ANGLE_RANGE_RWTH, getTomoAngleRangeRWTH(),	isTomoAngleRangeRWTH_defined()		);
		writeTag( t, TAG_DYN_RANGE_RWTH,		getDynRangeRWTH(),			isDynRangeRWTH_defined()			);
		writeTag( t, TAG_COL_END,				getCol_end(),				isCol_defined() 					);
		writeTag( t, TAG_COL_BEG,				getCol_beg(),				isCol_defined() 					);
		writeTag( t, TAG_ROW_END,				getRow_end(),				isRow_defined() 					);
		writeTag( t, TAG_ROW_BEG,				getRow_beg(),				isRow_defined() 					);

		if ( isCol_defined() )
		{
			if ( ( getCol_beg() < 0 ) || ( getCol_end() < 0 ) ||
				 ( (unsigned int )getCol_end() != (unsigned int )getCol_beg() + get_Dim_y() - 1 ) )
			{
				cerr << "WARNING. Inconsistency in col_end and col_beg." << endl;
				cerr << "col_end: " << getCol_end() << ", col_beg: " << getCol_beg() << endl;
			}
		}
		if ( isRow_defined() )
		{
			if ( ( getRow_beg() < 0 ) || ( getRow_end() < 0 ) ||
				 ( (unsigned int )getRow_end() != (unsigned int )getRow_beg() + get_Dim_x() - 1 ) )
			{
				cerr << "WARNING. Inconsistency in row_end and row_beg." << endl;
				cerr << "row_end: " << getRow_end() << ", row_beg: " << getRow_beg() << endl;
			}
		}
		
		for ( header_t::size_type tagcounter = 0; tagcounter < p_additionaltags.size(); tagcounter++ )
		{
			writeTag( t, p_additionaltags[tagcounter].tag, p_additionaltags[tagcounter].value );
		}

		t.push_back( "}" );

		// Fill in the number of lines of the header.
		writeTagIndex( t, TAG_HEADER_LINES_TUD, t.size(), 3 );

		// Now we are able to calculate the headersize.
		headersize = calculate_headersize( t, numCharsEndl, numSpacePadding );

		// Fill in the headersize value.
		writeTagIndex( t, TAG_HEADER_SIZE_TUD, headersize, 2 );

		// Fill in as many spaces as needed for a header size that is an interger
		// multiple of 2^10. Do not forget the closing curled bracket '}'.
		t[ t.size()-1 ] = string( numSpacePadding, ' ' ) + '}';

		// Start from i=1 because "{\n" has already been written to determine numCharsEndl.
		for ( text_t::size_type i=1; i<t.size(); i++ )
		{
			fprintf( f, "%s\n", t[i].c_str() );
		}

		switch ( getDataType() )
		{	
			case DT_DOUBLE:
			{
				if (BigEndianByteOrder())
					revByteOrder(mat, get_Dim_y()*get_Dim_x());
				if (flipByteOrder)
					revByteOrder(mat, get_Dim_y()*get_Dim_x());

				numBytesWritten = (unsigned long ) fwrite(mat, sizeof(double), get_Dim_x()*get_Dim_y(), f);

				if ( getVerbosity() )
					cout << numBytesWritten << " doubles written to " << filename << endl;

				if (BigEndianByteOrder())
					revByteOrder(mat, get_Dim_y()*get_Dim_x());
				if (flipByteOrder)
					revByteOrder(mat, get_Dim_y()*get_Dim_x());
			}
			break;

			case DT_FLOAT:
			{
				float *data = new float[get_Dim_x()*get_Dim_y()];

				for (unsigned int i = 0; i < get_Dim_y(); i++)
				{
					for (unsigned int j = 0; j < get_Dim_x(); j++)
					{
						data[((unsigned long )i*get_Dim_x()+(unsigned long )j)] = (float )mat[((unsigned long )i*get_Dim_x()+(unsigned long )j)];
						if (BigEndianByteOrder())
							revByteOrder(&(data[((unsigned long )i*get_Dim_x()+(unsigned long )j)]));
						if (flipByteOrder)
							revByteOrder(&(data[((unsigned long )i*get_Dim_x()+(unsigned long )j)]));
					}
				}		

				numBytesWritten = (unsigned long ) fwrite(data, sizeof(float), get_Dim_x()*get_Dim_y(), f);
				if ( getVerbosity() )
					cout << numBytesWritten << " floats written to " << filename << endl;

				delete []data;	
			}
			break;

			case DT_SHORT:
			{
				short *data;
				double range;
				double offset;
				double shift;

				data = new short[get_Dim_x()*get_Dim_y()];

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 - 1.0 );
					offset = get_ScalingMin();
					shift  = 65536.0 / 2.0;
				}
				else
				{
					range  = 1;
					offset = 0;
					shift  = 0;
				}

				for (unsigned int i = 0; i < get_Dim_y(); i++)
				{
					for (unsigned int j = 0; j < get_Dim_x(); j++)
					{
						index = (unsigned long )i * get_Dim_x() + (unsigned long )j;
						if ( range != 0 )
						{
							data[index] = (short )floor( (mat[index]-offset)/range-shift+epsilon );
						}
						else
						{
							data[index] = 0;
						}
					}
				}		

				if (flipByteOrder)
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());

				numBytesWritten = (unsigned long) fwrite(data, sizeof(short), get_Dim_x()*get_Dim_y(), f);
				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numBytesWritten << " shorts written to " << filename << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_UNSIGNED_SHORT:
			{
				unsigned short *data;
				double range;
				double offset;

				data = new unsigned short[get_Dim_x()*get_Dim_y()];

				if ( get_ScaleOut() )
				{
					range = ( get_ScalingMax() - get_ScalingMin() )/(double )65535;
					offset = get_ScalingMin();
				}
				else
				{
					range = 1;
					offset = 0;
				}

				for (unsigned int i = 0; i < get_Dim_y(); i++)
				{
					for (unsigned int j = 0; j < get_Dim_x(); j++)
					{
						index = (unsigned long )i * get_Dim_x() + (unsigned long )j;
						if ( range != 0 )
						{
							data[index] = (unsigned short )floor( (mat[index]-offset)/range+epsilon );
						}
						else
						{
							data[index] = 0;
						}
					}
				}		

				if (flipByteOrder)
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());

				numBytesWritten = (unsigned long) fwrite(data, sizeof(unsigned short), get_Dim_x()*get_Dim_y(), f);
				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numBytesWritten << " unsigned shorts written to " << filename << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_LONG:
			{
				long *data;
				double range;
				double offset;
				double shift;

				data = new long[get_Dim_x()*get_Dim_y()];

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 * 65536.0 - 1.0 );
					offset = get_ScalingMin();
					shift  = 65536.0 * 65536.0 / 2.0;
				}
				else
				{
					range  = 1;
					offset = 0;
					shift  = 0;
				}
				for (unsigned int i = 0; i < get_Dim_y(); i++)
				{
					for (unsigned int j = 0; j < get_Dim_x(); j++)
					{
						index = (unsigned long )i * get_Dim_x() + (unsigned long )j;
						if ( range != 0 )
						{
							data[index] = (long )floor( (mat[index]-offset)/range-shift+epsilon);
						}
						else
						{
							data[index] = 0;
						}
					}
				}		

				if (flipByteOrder)
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());

				numBytesWritten = (unsigned long) fwrite(data, sizeof(long), get_Dim_x()*get_Dim_y(), f);
				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numBytesWritten << " longs written to " << filename << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_UNSIGNED_LONG:
			{
				unsigned long *data;
				double range;
				double offset;

 				data = new unsigned long[get_Dim_x()*get_Dim_y()];

				if ( get_ScaleOut() )
				{
					range  = ( get_ScalingMax() - get_ScalingMin() ) / ( 65536.0 * 65536.0 - 1.0 );
					offset = get_ScalingMin();
				}
				else
				{
					range  = 1;
					offset = 0;
				}
				for (unsigned int i = 0; i < get_Dim_y(); i++)
				{
					for (unsigned int j = 0; j < get_Dim_x(); j++)
					{
						index = (unsigned long )i * get_Dim_x() + (unsigned long )j;
						if ( range != 0 )
						{
							data[index] = (unsigned long )floor( (mat[index]-offset)/range+epsilon);
						}
						else
						{
							data[index] = 0;
						}
					}
				}		

				if (flipByteOrder)
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_HIGH_BYTE_FIRST ) && ( !BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());
				if ( ( getByteOrder() == BO_LOW_BYTE_FIRST ) && ( BigEndianByteOrder() ) )
					revByteOrder(data, get_Dim_x()*get_Dim_y());

				numBytesWritten = (unsigned long) fwrite(data, sizeof(unsigned long), get_Dim_x()*get_Dim_y(), f);
				if ( getVerbosity() )
				{
					string un = ( get_ScaleOut() ? "" : "un" );
					cout << numBytesWritten << " unsigned longs written to " << filename << " (" << un << "scaled)" << endl;
				}

				delete []data;
			}
			break;

			case DT_NOT_SPECIFIED:
			{

			}
			break;

			case DT_FAMELESS:
			{

			}
			break;

			default:
			{
				// The user should never see this message.
				cerr << "ERROR! Bug in edf::get_datatype_str(). ";
				cerr <<	"Invalid data type '" << getDataType() << "'." << endl;
				setDataType( DT_NOT_SPECIFIED );
				retval = 1;
			}
			break;
		}

		fclose(f);
		return retval;
	}
*/
