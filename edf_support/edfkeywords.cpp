#include "edfkeywords.h"

#include <sstream>

#include "util.h"

using std::vector;
using std::string;
using std::ostringstream;
using std::ostream;
using matlibutil::equal;

namespace ns_key
{
	//------------------------------------------------------------------------------
	// Built-in keywords, i.e. keywords the values of which are interpreted by edf.
	//
	C_Keyword TAG_HEADER_ID(			"HeaderID" 								);
	C_Keyword TAG_TUD_HEADER_SIZE(		"TUD_HeaderSize"						);
	C_Keyword TAG_TUD_HEADER_LINES(		"TUD_HeaderLines"						);
	C_Keyword TAG_SIZE(					"Size"	 								);
	C_Keyword TAG_DIM(					"Dim_"	 								);
	C_Keyword TAG_BYTE_ORDER(			"ByteOrder"								);
	C_Keyword TAG_DATA_TYPE( 			"DataType" 								);
	C_Keyword TAG_REAL_FORMAT(			"RealFormat"							);
	C_Keyword TAG_TUD_SCALING_MIN(		"TUD_ScalingMin"	, "MinValRWTH"		);
	C_Keyword TAG_TUD_SCALING_MAX(		"TUD_ScalingMax"	, "MaxValRWTH"		);
	C_Keyword TAG_TUD_FILE_TYPE( 		"TUD_FileType"  	, "fileTypeRWTH"	);
	C_Keyword TAG_TUD_DATE_ORIG( 		"TUD_DateOrig" 							);
	C_Keyword TAG_TUD_DATE_PREV( 		"TUD_DatePrev" 							);
	C_Keyword TAG_TUD_DATE_CURR( 		"TUD_DateCurr" 							);
	C_Keyword TAG_TUD_FILENAME_ORIG( 	"TUD_FilenameOrig" 						);
	C_Keyword TAG_TUD_FILENAME_PREV( 	"TUD_FilenamePrev" 						);
	C_Keyword TAG_TUD_FILENAME_CURR( 	"TUD_FilenameCurr"						);
	C_Keyword TAG_IMAGE( 				"Image"									);
	C_Keyword TAG_VERSION_NUMBER(		"VersionNumber"		, "versionRWTH"		);

	//------------------------------------------------------------------------------
	// Recommended user-defined keywords.
	// User-defined keywords may be used by users, but they are not interpreted by
	// the edf engine. Recommended means, that the user should use those keywords
	// rather than keywords not in this list. If none of the recommended keywords
	// adequately describes the property to be stored in he header, the user may
	// define another one. If using a recommended user-defined keyword, the user
	// should not use the string literal itself, but the string constants defined
	// in the following list.
	//
	C_Keyword TAG_TUD_MATLIB_VERSION(	"TUD_MatlibVersion" 						);
	C_Keyword TAG_COMPRESSION(			"Compression"								);
	C_Keyword TAG_SR_CURRENT(			"SR_Current"								);
	C_Keyword TAG_TITLE( 				"Title" 									);
	C_Keyword TAG_PREFIX(				"prefix"									);
	C_Keyword TAG_SUFFIX(				"suffix"									);
	C_Keyword TAG_PRESET(				"preset"									);
	C_Keyword TAG_RUN(					"run"										);
	C_Keyword TAG_TUD_PATH(				"TUD_Path"									);
	C_Keyword TAG_TUD_FILENAME(			"TUD_Filename"								);
	C_Keyword TAG_TUD_BEG_X( 			"TUD_Beg_x"				, "col_beg"			);
	C_Keyword TAG_TUD_END_X( 			"TUD_End_x"				, "col_end"			);
	C_Keyword TAG_TUD_BEG_Y( 			"TUD_Beg_y"				, "row_beg"			);
	C_Keyword TAG_TUD_END_Y( 			"TUD_End_y"				, "row_end"			);
	C_Keyword TAG_TUD_BEG_Z( 			"TUD_Beg_z" 								);
	C_Keyword TAG_TUD_END_Z( 			"TUD_End_z" 								);
	C_Keyword TAG_TUD_STEP_X(			"TUD_Step_x"			, "xStepRWTH"		);
	C_Keyword TAG_TUD_STEP_Y(			"TUD_Step_y"			, "yStepRWTH"		);
	C_Keyword TAG_TUD_STEP_Z(			"TUD_Step_z"								);
	C_Keyword TAG_TUD_ZERO_PAD(			"TUD_ZeroPad"			, "zeroPadRWTH"		);
	C_Keyword TAG_TUD_ROT_AXIS(			"TUD_RotAxis"			, "rotAxisRWTH"		);
	C_Keyword TAG_TUD_ROT_AXIS_POS_X(	"TUD_RotAxisPosX"		, "rotAxisXPosRWTH"	);
	C_Keyword TAG_TUD_ROT_AXIS_POS_Y(	"TUD_RotAxisPosY"		, "rotAxisYPosRWTH"	);
	C_Keyword TAG_TUD_ROT_ANGLE( 		"TUD_RotAngle"			, "rotAngleRWTH"	);
	C_Keyword TAG_TUD_TOMO_ANGLE_RANGE(	"TUD_TomoAngleRange"	, "tomoAngleRangeRWTH"	);
	C_Keyword TAG_TUD_DYNAMIC_RANGE( 	"TUD_DynamicRange"		, "dynRangeRWTH"	);
	C_Keyword TAG_TUD_PROPOSAL_NUMBER(	"TUD_ProposalNumber"	, "ProposalNumber"	);
	C_Keyword TAG_TUD_SAMPLE_NAME(		"TUD_SampleName"		, "SampleName"		);
	C_Keyword TAG_TUD_DETECTOR_ALIAS( 	"TUD_DetectorAlias" 						);
	C_Keyword TAG_TUD_DETECTOR_NAME( 	"TUD_DetectorName" 		, "DetectorName"	);
	C_Keyword TAG_TUD_EXPOSURE_TIME( 	"TUD_ExposureTime"		, "count_time"		);
	C_Keyword TAG_TUD_BINNING_X( 		"TUD_Binning_x" 							);
	C_Keyword TAG_TUD_BINNING_Y( 		"TUD_Binning_y" 							);
	C_Keyword TAG_TUD_ORIENTATION(		"TUD_Orientation"							);
	C_Keyword TAG_TUD_TEMPERATURE(		"TUD_Temperature"							);
	C_Keyword TAG_TUD_ACCTR_DATE_CT( 	"TUD_AcctrDateCt"							);
	C_Keyword TAG_TUD_ACCTR_SETUP(		"TUD_AcctrSetup"							);
	C_Keyword TAG_TUD_ACCTR_SESSION( 	"TUD_AcctrSession"							);
	C_Keyword TAG_TUD_ACCTR_LOGFILE( 	"TUD_AcctrLogfile"		, "ScanFile"		);
	C_Keyword TAG_TUD_SCAN_NO(			"TUD_AcctrScanNo"		, "scan_no" 	, "ScanNumber"	);
	C_Keyword TAG_TUD_POINT_NO(			"TUD_AcctrPointNo"		, "point_no"	, "ScanPoint"	);
	C_Keyword TAG_TUD_MOTOR_NAMES(		"TUD_AcctrMotors"		, "motor_mne"		);
	C_Keyword TAG_TUD_MOTOR_POSITIONS(	"TUD_AcctrPositions"	, "motor_pos"		);
	C_Keyword TAG_TUD_COUNTER_NAMES( 	"TUD_AcctrCounters"		, "counter_mne"		);
	C_Keyword TAG_TUD_COUNTER_VALUES(	"TUD_AcctrCntValues"	, "counter_pos"		);

	//==========================================================================
	//================================================================ C_Keyword
	//
	//
	C_Keyword::C_Keyword()
	{
		set();
	}
	
	C_Keyword::C_Keyword( const string A )
	{
		set( A );
	}
	
	C_Keyword::C_Keyword( const string A, const string B )
	{
	   set( A, B );
	}
	
	C_Keyword::C_Keyword( const string A, const string B, const string C )
	{
		set( A, B, C );
	}

	void C_Keyword::set()
	{
		clear();
		push_back( "" );
	}
	
	void C_Keyword::set( const std::string A )
	{
		clear();
		push_back( A );
	}
	void C_Keyword::set( const std::string A, const std::string B )
	{
		clear();
		push_back( A );
		push_back( B );
	}
	void C_Keyword::set( const std::string A, const std::string B, const std::string C )
	{
		clear();
		push_back( A );
		push_back( B );
		push_back( C );
	}

	string C_Keyword::get_main() const
	{
		return *begin();
	}
	
	string C_Keyword::print() const
	{
		ostringstream osst;
		
		osst << "(";
		for ( vector<string>::const_iterator i = begin(); i != end(); i++ )
		{
			osst << "'" << *i << "'";
			if ( i+1 != end() )
			{
				osst << ", ";
			}
		}
		osst << ")";

		return osst.str();
	}

	bool operator==( const C_Keyword A, const C_Keyword B )
	{
		for ( vector<string>::const_iterator i = A.begin(); i != A.end(); i++ )
		{
			if ( *i == B )
			{
				return true;
			}
		}
		return false;
	}
	
	bool operator==( const C_Keyword A, const string b )
	{
		return b == A;
	}

	bool operator==( const string a, const C_Keyword B )
	{
		for ( vector<string>::const_iterator i = B.begin(); i != B.end(); i++ )
		{
			if ( equal( a, *i ) )
			{
				return true;
			}
		}
		return false;
	}

	ostream &operator<<( ostream &os, C_Keyword keyword )
	{
		return os << keyword.print();
	}
	
	


	C_AllKeywords::C_AllKeywords()
	{
		push_back( TAG_HEADER_ID			);	// Format: "EH:%6d:%6d:%6d", current, next, prvious.
		push_back( TAG_TUD_HEADER_SIZE		);	// Number of bytes of the header.
		push_back( TAG_TUD_HEADER_LINES 	);	// Number of header lines (inlcuding the lines with '{' and '}'.
		push_back( TAG_SIZE 				);	// Number of bytes of data (image data only, excluding header)
		push_back( TAG_DIM					);	// Number of pixels of directions 1, 2, 3, ...
		push_back( TAG_BYTE_ORDER			);	// Specifies the endianess LowByteFirst or HighByteFirst.
		push_back( TAG_DATA_TYPE			);	// Specifies datatype.
		push_back( TAG_REAL_FORMAT			);	// Only meaningful for float data types (IEEE, VAX, VMS, CONVEX, etc.)
		push_back( TAG_TUD_SCALING_MIN		);	// Min val for integer to double conversion.
		push_back( TAG_TUD_SCALING_MAX		);	// Min val for integer to double conversion.
		push_back( TAG_TUD_FILE_TYPE		);	// matrixdata, sinodata, etc.
		push_back( TAG_TUD_DATE_ORIG		);	// Date of image creation.
		push_back( TAG_TUD_DATE_PREV		);	// Date of previous modification.
		push_back( TAG_TUD_DATE_CURR		);	// Date of current modification.
		push_back( TAG_TUD_FILENAME_ORIG	);	// Name of the file when the image was created.
		push_back( TAG_TUD_FILENAME_PREV	);	// Name of the file the image had been loaded from.
		push_back( TAG_TUD_FILENAME_CURR	);	// Name of the file the image had been written to.
		push_back( TAG_IMAGE				);	// Image name (e.g. 1)
		push_back( TAG_VERSION_NUMBER		);	// Version of edf class writing the image.

		push_back( TAG_TUD_MATLIB_VERSION	);	// Matlib/Tomo version.
		push_back( TAG_COMPRESSION			);	// Compression method (jpg, DiffDataCompress, etc.)
		push_back( TAG_SR_CURRENT			);	// Synchrotron radiation current for normalization)
		push_back( TAG_TITLE				);	// ...to be further specified in future...
		push_back( TAG_PREFIX				);	//
		push_back( TAG_SUFFIX				);	//
		push_back( TAG_PRESET				);	//
		push_back( TAG_RUN					);	//
		push_back( TAG_TUD_PATH 			);	// Directory path of the image when stored for the 1st time.
		push_back( TAG_TUD_FILENAME 		);	// Filename of the image when stored for the 1st time.
		push_back( TAG_TUD_BEG_X			);	// First x-coordinate of a ROI (including this value).
		push_back( TAG_TUD_END_X			);	// Last x-coordinate of a ROI (excuding this value).
		push_back( TAG_TUD_BEG_Y			);	// First y-coordinate of a ROI (including this value).
		push_back( TAG_TUD_END_Y			);	// Last y-coordinate of a ROI (excuding this value).
		push_back( TAG_TUD_BEG_Z			);	// First z-coordinate of a ROI (including this value).
		push_back( TAG_TUD_END_Z			);	// Last z-coordinate of a ROI (excuding this value).
		push_back( TAG_TUD_STEP_X			);	// Pixelsize in x-direction.
		push_back( TAG_TUD_STEP_Y			);	// Pixelsize in y-direction.
		push_back( TAG_TUD_STEP_Z			);	// Pixelsize in z-direction.
		push_back( TAG_TUD_ZERO_PAD 		);	// Faktor for zero padding size.
		push_back( TAG_TUD_ROT_AXIS 		);	// Rotation axis for tomographic rec.
		push_back( TAG_TUD_ROT_AXIS_POS_X	);	// Rotation axis x-pos in image.
		push_back( TAG_TUD_ROT_AXIS_POS_Y	);	// Rotation axis y-pos in image.
		push_back( TAG_TUD_ROT_ANGLE		);	// Integer rotation angle for tomographic rec.
		push_back( TAG_TUD_TOMO_ANGLE_RANGE );	// Angle range for tomography in multiples of 180.
		push_back( TAG_TUD_DYNAMIC_RANGE	);	// Dynamic range of image data.
		push_back( TAG_TUD_PROPOSAL_NUMBER	);	// ESRF proposal number.
		push_back( TAG_TUD_SAMPLE_NAME		);	// Sample name.
		push_back( TAG_TUD_DETECTOR_ALIAS	);	// Name of the detector the image was exposed with.
		push_back( TAG_TUD_DETECTOR_NAME	);	// Name of the detector the image was exposed with.
		push_back( TAG_TUD_EXPOSURE_TIME	);	// Exposure time in seconds.
		push_back( TAG_TUD_BINNING_X		);	// Detector binning in x-direction.
		push_back( TAG_TUD_BINNING_Y		);	// Detector binning in y-direction.
		push_back( TAG_TUD_ORIENTATION		);	// Detector orientation (normal, rot90, rot180, rot270, mirrorUpDown, mirrorLeftRight)
		push_back( TAG_TUD_TEMPERATURE		);	// Temperature of ccd chip.
		push_back( TAG_TUD_ACCTR_DATE_CT	);	// Date and time of (beginning of) detector exposure.
		push_back( TAG_TUD_ACCTR_SETUP		);	// Acctr setup name (in spec this would be the session name).
		push_back( TAG_TUD_ACCTR_SESSION	);	// Acctr session name (within a setup).
		push_back( TAG_TUD_ACCTR_LOGFILE	);	// Acctr logfile name.
		push_back( TAG_TUD_SCAN_NO			);	// Acctr scan number (within a session).
		push_back( TAG_TUD_POINT_NO 		);	// Acctr scan point number (within a scan).
		push_back( TAG_TUD_MOTOR_NAMES		);	// Motor names defined in acctr setup (space separated list).
		push_back( TAG_TUD_MOTOR_POSITIONS	);	// Positions of all motors (space separated list).
		push_back( TAG_TUD_COUNTER_NAMES	);	// Counter names defined in acctr setup (space separated list).
		push_back( TAG_TUD_COUNTER_VALUES	);	// Values of all motors (space separated list).
	}
}
