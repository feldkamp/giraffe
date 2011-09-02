#ifndef edfkeywords_h_flag
#define edfkeywords_h_flag

#include <vector>
#include <string>

namespace ns_key
{
	//================================================================ C_Keyword
	//
	//
	class C_Keyword : public std::vector<std::string>
	{
		public:
			C_Keyword();
			C_Keyword( const std::string A );
			C_Keyword( const std::string A, const std::string B );
			C_Keyword( const std::string A, const std::string B, const std::string C );
			
			void set();
			void set( const std::string A );
			void set( const std::string A, const std::string B );
			void set( const std::string A, const std::string B, const std::string C );
			
			std::string get_main() const;
			
			std::string print() const;
	};
	bool operator==( const C_Keyword A, const C_Keyword B );
	bool operator==( const C_Keyword A, const std::string B );
	bool operator==( const std::string a, const C_Keyword B );
	std::ostream &operator<<( std::ostream &os, C_Keyword keyword );



	//------------------------------------------------------------------------------
	// Built-in keywords, i.e. keywords the values of which are interpreted by edf.
	//
	extern C_Keyword TAG_HEADER_ID;
	extern C_Keyword TAG_TUD_HEADER_SIZE;
	extern C_Keyword TAG_TUD_HEADER_LINES;
	extern C_Keyword TAG_SIZE;
	extern C_Keyword TAG_DIM;
	extern C_Keyword TAG_BYTE_ORDER;
	extern C_Keyword TAG_DATA_TYPE;
	extern C_Keyword TAG_REAL_FORMAT;
	extern C_Keyword TAG_TUD_SCALING_MIN;
	extern C_Keyword TAG_TUD_SCALING_MAX;
	extern C_Keyword TAG_TUD_FILE_TYPE;
	extern C_Keyword TAG_TUD_DATE_ORIG;
	extern C_Keyword TAG_TUD_DATE_PREV;
	extern C_Keyword TAG_TUD_DATE_CURR;
	extern C_Keyword TAG_TUD_FILENAME_ORIG;
	extern C_Keyword TAG_TUD_FILENAME_PREV;
	extern C_Keyword TAG_TUD_FILENAME_CURR;
	extern C_Keyword TAG_IMAGE;
	extern C_Keyword TAG_VERSION_NUMBER;

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
	extern C_Keyword TAG_TUD_MATLIB_VERSION;
	extern C_Keyword TAG_COMPRESSION;
	extern C_Keyword TAG_SR_CURRENT;
	extern C_Keyword TAG_TITLE;
	extern C_Keyword TAG_PREFIX;
	extern C_Keyword TAG_SUFFIX;
	extern C_Keyword TAG_PRESET;
	extern C_Keyword TAG_RUN;
	extern C_Keyword TAG_TUD_PATH;
	extern C_Keyword TAG_TUD_FILENAME;
	extern C_Keyword TAG_TUD_BEG_X;
	extern C_Keyword TAG_TUD_END_X;
	extern C_Keyword TAG_TUD_BEG_Y;
	extern C_Keyword TAG_TUD_END_Y;
	extern C_Keyword TAG_TUD_BEG_Z;
	extern C_Keyword TAG_TUD_END_Z;
	extern C_Keyword TAG_TUD_STEP_X;
	extern C_Keyword TAG_TUD_STEP_Y;
	extern C_Keyword TAG_TUD_STEP_Z;
	extern C_Keyword TAG_TUD_ZERO_PAD;
	extern C_Keyword TAG_TUD_ROT_AXIS;
	extern C_Keyword TAG_TUD_ROT_AXIS_POS_X;
	extern C_Keyword TAG_TUD_ROT_AXIS_POS_Y;
	extern C_Keyword TAG_TUD_ROT_ANGLE;
	extern C_Keyword TAG_TUD_TOMO_ANGLE_RANGE;
	extern C_Keyword TAG_TUD_DYNAMIC_RANGE;
	extern C_Keyword TAG_TUD_PROPOSAL_NUMBER;
	extern C_Keyword TAG_TUD_SAMPLE_NAME;
	extern C_Keyword TAG_TUD_DETECTOR_ALIAS;
	extern C_Keyword TAG_TUD_DETECTOR_NAME;
	extern C_Keyword TAG_TUD_EXPOSURE_TIME;
	extern C_Keyword TAG_TUD_BINNING_X;
	extern C_Keyword TAG_TUD_BINNING_Y;
	extern C_Keyword TAG_TUD_ORIENTATION;
	extern C_Keyword TAG_TUD_TEMPERATURE;
	extern C_Keyword TAG_TUD_ACCTR_DATE_CT;
	extern C_Keyword TAG_TUD_ACCTR_SETUP;
	extern C_Keyword TAG_TUD_ACCTR_SESSION;
	extern C_Keyword TAG_TUD_ACCTR_LOGFILE;
	extern C_Keyword TAG_TUD_SCAN_NO;
	extern C_Keyword TAG_TUD_POINT_NO;
	extern C_Keyword TAG_TUD_MOTOR_NAMES;
	extern C_Keyword TAG_TUD_MOTOR_POSITIONS;
	extern C_Keyword TAG_TUD_COUNTER_NAMES;
	extern C_Keyword TAG_TUD_COUNTER_VALUES;

	class C_AllKeywords: public std::vector<C_Keyword>
	{
		public:
			C_AllKeywords();
	};
}


#endif
