#ifndef ANALYZE_H
#define ANALYZE_H
//
//  analyze.h
//  xcca_commandline
//
//  Created by Feldkamp on 7/21/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "arrayclasses.h"
#include "arraydataIO.h"

class Analyzer {
public:	
	Analyzer();
	~Analyzer();

	int processFiles( vector<string> files, double cenX, double cenY, int num_phi, int num_q);
	int processFiles( vector<string> files, double cenX, double cenY, int num_phi, int num_q, 
						double start_q, double stop_q, int LUTx, int LUTy );
						
	void setBackground( array2D *back );
	array2D *background();
	void setBackgroundWeight( double weight );
	double backgroundWeight();
	void setMask( array2D *newmask );
	array2D *mask();
	void setOutputDirectory( string outdir );
	string outputDirectory();
	void setAlg( int alg );
	int alg();
	
	
	void setFlagSubtractBackground( bool flag );
	bool flagSubtractBackground();
	void setFlagSingleCorrelationOutput( bool flag );
	bool flagSingleCorrelationOutput();
	void setFlagUseMask( bool flag );
	bool flagUseMask();
	
private:
	arraydataIO *io;
	array2D *p_back;
	array2D *p_mask;
	double p_back_weight;
	string p_out_dir;
	int p_alg;
	bool p_flag_subtract_background;				// subtract background?
	bool p_flag_single_correlation_output;		// write every single correlation?
	bool p_flag_use_mask;
};

#endif
