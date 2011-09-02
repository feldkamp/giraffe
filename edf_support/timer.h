#ifndef timer_h_flag
#define timer_h_flag

#ifdef WIN32
	#pragma once
	#include <winsock2.h>			// by JF: winsock2 is not really needed in timer, but if it is at all used somewhere in the project,
									// it absolutely MUST be included before windows.h, this is a safety precaution
	#include <windows.h>
	
	typedef unsigned int useconds_t;	//useconds_t is not known on windows machines
#endif

#include <string>

namespace ns_timer
{
	enum unit_t
	{
		none=0,
		usec=1,
		msec=2,
		sec=3
	};

	enum date_format_t
	{
		DATE_FORMAT_FILENAME,
		DATE_FORMAT_GERMAN,
		DATE_FORMAT_AMERICAN,
		DATE_FORMAT_AMERICAN_2
	};

	enum date_entries_t
	{
		DATE_ENTRIES_YEAR		= 1L << 1,
		DATE_ENTRIES_MONTH		= 1L << 2,
		DATE_ENTRIES_DAY		= 1L << 3,
		DATE_ENTRIES_HOUR		= 1L << 4,
		DATE_ENTRIES_MINUTES	= 1L << 5,
		DATE_ENTRIES_SECONDS	= 1L << 6,
		DATE_ENTRIES_MSEC		= 1L << 7,
		DATE_ENTRIES_USEC		= 1L << 8,
		
		DATE_ENTRIES_ALL = DATE_ENTRIES_YEAR |
						   DATE_ENTRIES_MONTH |
						   DATE_ENTRIES_DAY |
						   DATE_ENTRIES_HOUR |
						   DATE_ENTRIES_MINUTES |
						   DATE_ENTRIES_SECONDS |
						   DATE_ENTRIES_MSEC |
						   DATE_ENTRIES_USEC,

		DATE_ENTRIES_DEFAULT = DATE_ENTRIES_ALL & (~DATE_ENTRIES_MSEC) & (~DATE_ENTRIES_USEC)
	};

	//=================================================================== C_Date
	//
	class C_Date
	{
		public:
			C_Date();
			C_Date( int year, int month, int day,
					int hour, int minutes, int seconds,
					int msec=0, int usec=0 );
			
			void init();
			
			template<class T>
			void set( T year, T month, T day,
					  T hour, T minutes, T seconds,
					  T msec=0, T usec=0 )
			{
				set_int( (int )year, (int )month, (int )day,
						 (int )hour, (int )minutes, (int )seconds,
						 (int )msec, (int )usec );
			};
			void set_int( int year, int month, int day,
						  int hour, int minutes, int seconds,
						  int msec=0, int usec=0 );
			void set( const C_Date &date );
			void set_year( int year );
			void set_month( int month );
			void set_day( int day );
			void set_hour( int hour );
			void set_minutes( int minutes );
			void set_seconds( int seconds );
			void set_msec( int msec );
			void set_usec( int usec );
			void set_date_now();
			void set_date_format( date_format_t date_format );
			void set_date_entries( date_entries_t date_entries );

			int get_year() const;
			int get_month() const;
			int get_day() const;
			int get_hour() const;
			int get_minutes() const;
			int get_seconds() const;
			int get_msec() const;
			int get_usec() const;
			date_format_t get_date_format() const;
			date_entries_t get_date_entries() const;
			C_Date get_date_now() const;
			C_Date get_date_now( date_format_t date_format ) const;
			std::string toMonthName( int number ) const;
			std::string print() const;
			
		private:
			int p_year;
			int p_month;
			int p_day;
			int p_hour;
			int p_minutes;
			int p_seconds;
			int p_msec;
			int p_usec;
			date_format_t p_date_format;
			date_entries_t p_date_entries;
	};
	
	std::ostream &operator<<( std::ostream &os, const C_Date &date );


	//================================================================== C_Timer
	//
	class C_Timer
	{
		public:
			C_Timer();
			
			double time_now( unit_t unit ) const;
			void set_timer( double time, unit_t unit );
			double get_timer( unit_t unit ) const;
			double get_remaining_time( unit_t unit ) const;
			double get_elapsed_time( unit_t unit ) const;
			bool timeout() const;
			bool activated() const;
			void reset();
			
		private:
			void init_time();
			double convert( double timeval, unit_t unit1, unit_t unit2 ) const;
			unit_t get_unit() const;
			#ifndef WIN32					
				double p_time_sec_0;				// Time in seconds when the timer object was created.
			#else
				LARGE_INTEGER p_time_largeint_0;	// Time in 100ns when the timer object was created.
			#endif
			double p_time_start;
			double p_time_delta;
			bool   p_activated;
			unit_t p_unit;
	};

	#ifdef WIN32 
		//----------------------------------------------------------------------
		// On windows, the time is stored in large integers (64-bit).
		// They represent 100 nano-seconds, but the actual accuracy
		// is limited to milliseconds..

		//------------------------------------------------- get_windowstimeofday
		//
		// Returns the system time as a large (64-bit) integer.
		//
		LARGE_INTEGER get_windowstimeofday();
		
		//--------------------------------------------------- get_windowstimeStr
		//
		// Return the human-readable system time on windows machines.
		//
		std::string get_windowstimeStr();		

		//---------------------------------------------------------------- alarm
		//
		// Define the alarm() function on windows machines.
		// (On unix included with unistd.h.)
		// Sleeps for 'seconds' and issues a signal SIGINT after that.
		//
		unsigned int alarm(unsigned int seconds);
	#endif

	//--------------------------------------------------------------- usleep
	//
	// Implements the usleep() function on windows machines.
	// Sleeps for 'time_usec' micro seconds.
	//
	int usleep( useconds_t time_usec );
	
} // namespace

#endif	// timer_h_flag
