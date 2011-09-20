//
//  main.cpp
//  TestSuite
//
//  Created by Feldkamp on 8/23/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include "tester.h"



void usage(){
	cout << "====================usage============================" << endl;
	cout << " -t<num>        : test different features" << endl;
	cout << "                      num=1 : testCrossCorrelator" << endl;	
	cout << "                      num=2 : testArrayClasses" << endl;	
	cout << "                      num=3 : testFourierTrafo" << endl;	
	cout << "                      num=4 : testIO" << endl;	
	cout << "                      num=5 : testDataTypes" << endl;	
	cout << "                      num=0 : test all" << endl;	
	cout << "=====================================================" << endl;
}

string argToString( char * const argument ){
	//fill argument after the -<option letter> into a string
	std::ostringstream osst;
	int j = 2;
	while ( argument[j] != '\0' ){
		osst << argument[j];
		j++;
	}
	string retstring = osst.str();
	return retstring;
}

int main (int argc, char * const argv[]) {

    cout << "Hello, World!\n";

	if (argc < 2){
		usage();
		return 1;
	} else {
		for(int i = 1; i < argc; i++){						// check if all options are valid first
			if (argv[i][0] != '-')
			{
				cout << "'" << argv[i] << "' does not start with a dash: not a valid option. Exiting." << endl;
				usage();
				return 1;
			}
		}

		for(int i = 1; i < argc; i++){						//if all options are valid, proceed and evaluate
			
			if (argv[i][1] == 't') {
				cout << "---- testing option '" << argv[i][2] << "' ----" << endl;
				//---------------------------------------------------
				//run various tests
				//---------------------------------------------------
				string base = "/Users/feldkamp/Desktop/";
				cout << "output directory '" << base << "'" << endl;
				Tester *t = new Tester();
				t->setBase(base);
				
				switch(argv[i][2]){
					case '1':
						t->testCrossCorrelator( 1 );					
						break;
					case '2':
						t->testArrayClasses();
						break;
					case '3':
						t->testFourierTrafo();
						break;
					case '4':
						t->testIO();
						break;
					case '5':
						t->testDataTypes();					
						break;
					case '6':
						t->testArraySpeed();
						break;
					case '0':							// fall through to default
					default:
						t->testCrossCorrelator( 1 );
						t->testArrayClasses();
						t->testFourierTrafo();
						t->testIO();
						t->testDataTypes();
						break;
				}//end switch
				cout << "---- testing done ----" << endl;
				delete t;
				return 0;
			} else {
				cout << "-" << argv[i][1] << " is not a valid option." << endl;
				usage();
				return 2;
			}
		}//end for i
	}//end if
	
	
	return 0;
}