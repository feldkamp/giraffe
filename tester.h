//
//  tester.h
//  TestSuite
//
//  Created by Feldkamp on 8/23/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//


#include <string>


class Tester {
	public:
		Tester();
		Tester( std::string base );
		~Tester();
		
		void setBase( std::string base );
		std::string base();

		int testArrayClasses();												// case 1
		int testCrossCorrelator( int alg = 1, int testpattern = 1 );		// case 2
		int testFourierTrafo();												// case 3
		int testIO(int mode=0);												// case 4 (mode 0: all formats, 1:edf, 2:hdf5, 3:tiff, 4:ascii)
		int testDataTypes();												// case 5
		int testArraySpeed();												// case 6
		
		//sub-tests for testArraySpeed
		void subtestBoost();
		void subtestArraydata(unsigned int size);
		void subtestVector(unsigned int size);
		void subtest2D_a(unsigned int dim1, unsigned int dim2);
		void subtest2D_b(unsigned int dim1, unsigned int dim2);
			
	private:
		std::string p_base;
};

