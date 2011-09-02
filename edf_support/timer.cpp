#include "timer.h"

#ifdef WIN32
	#pragma once
	#include <ctime>
	#include <cmath>
#else
	#include <sys/time.h>
#endif

#include <iostream>
#include <sstream>
#include <csignal>				//needed to 'raise(signal)' in 'alarm()'

#include "util.h"

using std::cout;
using std::endl;
using std::ostringstream;
using matlibutil::appendNumber;


namespace ns_timer
{
	//=================================================================== C_Date
	//
	//
	C_Date::C_Date()
	{
		init();
	}
	
	C_Date::C_Date( int year, int month, int day,
					int hour, int minutes, int seconds,
					int msec, int usec )
	{
		init();
		set( year, month, day, hour, minutes, seconds, msec, usec );
	}
	
	void C_Date::init()
	{
		set( 0, 0, 0, 0, 0, 0, 0, 0 );
		p_date_format = DATE_FORMAT_GERMAN;
		p_date_entries = DATE_ENTRIES_DEFAULT;
	}

	void C_Date::set_int( int year, int month, int day,
						  int hour, int minutes, int seconds,
						  int msec, int usec )
	{
		set_year( year );
		set_month( month );
		set_day( day );
		set_hour( hour );
		set_minutes( minutes );
		set_seconds( seconds );
		set_msec( msec );
		set_usec( usec );
	}
	
	void C_Date::set( const C_Date &date )
	{
		*this = date;
	}
	
	void C_Date::set_year( int year )
	{
		p_year = year;
	}
	
	void C_Date::set_month( int month )
	{
		p_month = month;
	}
	
	void C_Date::set_day( int day )
	{
		p_day = day;
	}
	
	void C_Date::set_hour( int hour )
	{
		p_hour = hour;
	}
	
	void C_Date::set_minutes( int minutes )
	{
		p_minutes = minutes;
	}
	
	void C_Date::set_seconds( int seconds )
	{
		p_seconds = seconds;
	}
	
	void C_Date::set_msec( int msec )
	{
		p_msec = msec;
	}
	
	void C_Date::set_usec( int usec )
	{
		p_usec = usec;
	}
	
	void C_Date::set_date_now()
	{
		set( get_date_now() );
	}

	void C_Date::set_date_format( date_format_t date_format )
	{
		p_date_format = date_format;
	}

	void C_Date::set_date_entries( date_entries_t date_entries )
	{
		p_date_entries = date_entries;
	}

	int C_Date::get_year() const
	{
		return p_year;
	}
	
	int C_Date::get_month() const
	{
		return p_month;
	}
	
	int C_Date::get_day() const
	{
		return p_day;
	}
	
	int C_Date::get_hour() const
	{
		return p_hour;
	}
	
	int C_Date::get_minutes() const
	{
		return p_minutes;
	}
	
	int C_Date::get_seconds() const
	{
		return p_seconds;
	}

	int C_Date::get_msec() const
	{
		return p_msec;
	}

	int C_Date::get_usec() const
	{
		return p_usec;
	}

	date_format_t C_Date::get_date_format() const
	{
		return p_date_format;
	}

	date_entries_t C_Date::get_date_entries() const
	{
		return p_date_entries;
	}

	C_Date C_Date::get_date_now( date_format_t date_format ) const
	{
		C_Date date;

		date.set_date_format( date_format );
		return date.get_date_now();
	}

	C_Date C_Date::get_date_now() const
	{
		C_Date date;

		date = *this;

		struct tm *mydate;	// formated time info (dd,mmm,yy,h,m,s)
		time_t timestamp;	// seconds since January 1, 1970

		timestamp = time(0);		// for date	
		
		mydate = localtime( &timestamp );		//make corrections from coordinated universal time (UTC) 

		date.set_year( mydate->tm_year + 1900 );
		date.set_month( mydate->tm_mon + 1 );
		date.set_day( mydate->tm_mday );
		date.set_hour( mydate->tm_hour );
		date.set_minutes( mydate->tm_min );
		date.set_seconds( mydate->tm_sec );
			
		#ifndef WIN32
			struct timeval tod;					// for usec
			gettimeofday( &tod, 0 );			// for usec
			date.set_msec( tod.tv_usec / 1000 );
			date.set_usec( tod.tv_usec % 1000 );
		#else
			date.set_msec( 0 ); // <--------- Still has to be implemented.
			date.set_usec( 0 ); // <--------- Still has to be implemented.
		#endif
		
		//old WINDOWS code
/*			SYSTEMTIME st;			

			GetSystemTime(&st);			

			date.set_year( st.wYear );
			date.set_month( st.wMonth );
			date.set_day( st.wDay );
			date.set_hour( st.wHour );
			date.set_minutes( st.wMinute );
			date.set_seconds( st.wSecond );
			date.set_msec( 0 ); // <--------- Still has to be implemented.
			date.set_usec( 0 ); // <--------- Still has to be implemented.
*/
		
		return date;
	}
	
	std::string C_Date::toMonthName( int number ) const
	{
		switch ( number )
		{
			case  1: return "Jan";
			case  2: return "Feb";
			case  3: return "Mar";
			case  4: return "Apr";
			case  5: return "Mai";
			case  6: return "Jun";
			case  7: return "Jul";
			case  8: return "Aug";
			case  9: return "Sep";
			case 10: return "Oct";
			case 11: return "Nov";
			case 12: return "Dec";
			default:
				cout << "ERROR. In C_Date::toMonthName() - "
					<< "Invalid month number " << number << "." << endl;
				return "---";
		}
	}

	std::string C_Date::print() const
	{
		ostringstream osst;

		switch ( get_date_format() )
		{
			case DATE_FORMAT_FILENAME:
			{
				osst << appendNumber( "", get_year(), 4 ) << "-"
					 << appendNumber( "", get_month(), 2 ) << "-"
					 << appendNumber( "", get_day(), 2 ) << "_";
				osst << appendNumber( "", get_hour(), 2 ) << "h"
					 << appendNumber( "", get_minutes(), 2 ) << "m"
					 << appendNumber( "", get_seconds(), 2 ) << "s";
					 
				if ( get_date_entries() & DATE_ENTRIES_MSEC )
				{
					osst << "." << appendNumber( "", get_hour(), 3 );
				}
			} break;
			
			case DATE_FORMAT_GERMAN:
			{
				osst << appendNumber( "", get_day(), 2 ) << "."
					 << appendNumber( "", get_month(), 2 ) << "."
					 << appendNumber( "", get_year(), 4 ) << ", ";
				osst << appendNumber( "", get_hour(), 2 ) << ":"
					 << appendNumber( "", get_minutes(), 2 ) << ":"
					 << appendNumber( "", get_seconds(), 2 );
					 
				if ( get_date_entries() & DATE_ENTRIES_MSEC )
				{
					osst << "." << appendNumber( "", get_hour(), 3 );
				}
			} break;
			
			case DATE_FORMAT_AMERICAN:
			{
				osst << appendNumber( "", get_month(), 2 ) << "/"
					 << appendNumber( "", get_day(), 2 ) << "/"
					 << appendNumber( "", get_year(), 4 ) << ", ";
				osst << appendNumber( "", get_hour(), 2 ) << ":"
					 << appendNumber( "", get_minutes(), 2 ) << ":"
					 << appendNumber( "", get_seconds(), 2 );

				if ( get_date_entries() & DATE_ENTRIES_MSEC )
				{
					osst << "." << appendNumber( "", get_hour(), 3 );
				}
			} break;
			
			case DATE_FORMAT_AMERICAN_2:
			{
				osst << toMonthName( get_month() ) << " "
					 << appendNumber( "", get_day(), 2 ) << " "
					 << appendNumber( "", get_year(), 4 ) << ", ";
				osst << appendNumber( "", get_hour(), 2 ) << ":"
					 << appendNumber( "", get_minutes(), 2 ) << ":"
					 << appendNumber( "", get_seconds(), 2 );

				if ( get_date_entries() & DATE_ENTRIES_MSEC )
				{
					osst << "." << appendNumber( "", get_hour(), 3 );
				}
			} break;
			
			default:
			{
				cout << "WARNING. Invalid date format '" << get_date_format() << "'." << endl;
			} break;
		}
		return osst.str();
	}

	std::ostream &operator<<( std::ostream &os, const class C_Date &date )
	{
		return os << date.print();
	}



	//================================================================== C_Timer
	//
	//
	C_Timer::C_Timer()
	{
		p_unit = usec; // unit for internal time values
		
		// Ôinit_time()Õ has been integrated in Ôreset()Õ. (JP, Dec08)
		//init_time();
		reset();
	}

	double C_Timer::time_now( unit_t unit ) const
	{
		double time;
		
#ifndef WIN32
		timeval time_now;
		gettimeofday( &time_now, 0 );
		
		switch ( unit )
		{
			case usec: time = 1000*1000*(time_now.tv_sec - p_time_sec_0) + time_now.tv_usec; break;
			case msec: time = 1000*(time_now.tv_sec - p_time_sec_0) + time_now.tv_usec/1000.; break;
			case sec:  time = (time_now.tv_sec - p_time_sec_0) + time_now.tv_usec/(1000.*1000.); break;
			case none: time = 0; break;
			default:   time = 0;
		}
#else	//-----------------------------------------------------------------WINDOWS CASE
		LARGE_INTEGER time_largeint_now = get_windowstimeofday();
		double time_usec_passed = (double) (time_largeint_now.QuadPart - p_time_largeint_0.QuadPart)/10.;
		//the factor 10 is needed to convert to mu-seconds
		//ex. for a LARGE_INTEGER: 128304466173430000
		//would be                 12830446617343000  usecs
		//or                       12830446617343     msecs
		
		switch ( unit )
		{
			case usec: time = time_usec_passed; break;
			case msec: time = time_usec_passed/1000.; break;
			case sec:  time = time_usec_passed/(1000.*1000.); break;
			case none: time = 0; break;
			default:   time = 0;
		}
#endif
		return time;
	}

	void C_Timer::set_timer( double time, unit_t unit )
	{
		reset();
		p_time_start = time_now( get_unit() );
		p_time_delta = convert( time, unit, get_unit() );
		p_activated = true;
	}

	double C_Timer::get_timer( unit_t unit ) const
	{
		return convert( p_time_delta, get_unit(), unit );
	}
	
	double C_Timer::get_remaining_time( unit_t unit ) const
	{		
		return convert( p_time_start + p_time_delta - time_now( get_unit() ), get_unit(), unit );
	}

	double C_Timer::get_elapsed_time( unit_t unit ) const
	{
		return convert( time_now( get_unit() ) - p_time_start, get_unit(), unit );
	}

	bool C_Timer::timeout() const
	{
		return (( get_remaining_time( get_unit() ) <= 0 ) && activated());
	}
	
	void C_Timer::reset()
	{
		// The following line is new and has not been tested, yet. (JP, Dec08)
		init_time();
		
		p_time_start = time_now( get_unit() );
		p_time_delta = 0;
		p_activated = false;
	}
	

	void C_Timer::init_time()
	{
#ifndef WIN32
		timeval time_now;
		gettimeofday( &time_now, 0 );
		// The following change has not been tested, yet. (JP, Dec08)
		//p_time_sec_0 = time_now.tv_sec;
		p_time_sec_0 = time_now.tv_sec + time_now.tv_usec/1000./1000.;
#else
		p_time_largeint_0 = get_windowstimeofday(); 
#endif	
	}

	double C_Timer::convert( double timeval, unit_t unit1, unit_t unit2 ) const
	{
		double result;
		
		if ( unit1 == unit2 )
		{
			result = timeval;
		}
		else if ( ( unit1 == usec ) && ( unit2 == msec ) ||
				  ( unit1 == msec ) && ( unit2 == sec ) )
		{
			result = timeval / 1000;
		}
		else if ( ( unit2 == usec ) && ( unit1 == msec ) ||
				  ( unit2 == msec ) && ( unit1 == sec ) )
		{
			result = timeval * 1000;
		}
		else if ( ( unit1 == usec ) && ( unit2 == sec ) )
		{
			result = timeval / (1000*1000);
		}
		else if ( ( unit2 == usec ) && ( unit1 == sec ) )
		{		
			result = timeval * (1000*1000);
		}
		else
		{
			result = 0;
		}
	
		return result;
	}
	
	bool C_Timer::activated() const
	{
		return p_activated;
	}
	
	unit_t C_Timer::get_unit() const
	{
		return p_unit;
	}
	
//--------------------------------------------------------------end of C_Timer




#ifdef WIN32
	// define usleep() on windows machines so that it can be used consistently for WIN32 and linux 
	// by specifying the namespace ns_timer whenever calling usleep;
	int usleep( useconds_t time_usec ) 
	{
		Sleep( (DWORD )ceil( (double )time_usec/1000.0 ) );
		return 0;
	}
	
	LARGE_INTEGER get_windowstimeofday()
	{
		SYSTEMTIME st;			
		FILETIME ft;
		LARGE_INTEGER li;		//64-bit signed integer 2^64 = 1.84E+19
								//this counts in tenths of usecs, so it can hold 2^64/2 'ticks'
								//2^63 / (10*1000*1000*60*60*24*365) = 29247 years

		GetSystemTime(&st);					//get SYSTEM
		SystemTimeToFileTime(&st, &ft);		//convert to FILETIME
    
		li.LowPart = ft.dwLowDateTime;		//write to LARGE_INTEGER
		li.HighPart = ft.dwHighDateTime;
		//cout << "li.QuadPart: " << li.QuadPart << endl;
		return li;
	}

	//display and return the current system time on windows machines
	std::string get_windowstimeStr()
	{	
		std::string outputstring = "";
		SYSTEMTIME st;			
		GetSystemTime(&st);			

		outputstring = "" + matlibutil::toString(st.wYear) + "/" + matlibutil::toString(st.wMonth) + "/" + matlibutil::toString(st.wDay)
		+ ", " + matlibutil::toString(st.wHour) + ":" + matlibutil::toString(st.wMinute) + ":" + matlibutil::toString(st.wSecond) 
		+ ", " + matlibutil::toString(st.wMilliseconds) + "ms";
		
		return outputstring;
	}
	
	//sleep for 'seconds' and issue a signal SIGINT (ctrl-c) after that
	unsigned int alarm(unsigned int seconds)
	{
		ns_timer::usleep( seconds*1000*1000 );
		return (unsigned int) raise( SIGINT );
	}
	
#else	// linux and macintosh have the 'usleep()' routine built-in
		// get_WINtime... routines can be used only on windows machines
	int usleep( useconds_t time_usec )
	{
		// Important note: Do not omit the double colon to call the standard usleep()!
		// Infinite recursive calling of functions may be annoying...
		::usleep(time_usec);
		return 0;
	}
#endif

} // namespace ns_timer
