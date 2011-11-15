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
#include <vector>

#include "arrayclasses.h"
#include "arraydataIO.h"

class Analyzer {
public:	
	Analyzer();
	~Analyzer();

	int processFiles( std::vector<std::string> files, int num_phi, int num_q);
	int processFiles( std::vector<std::string> files, int num_phi, int num_q, double start_q, double stop_q, int LUTx, int LUTy, 
							std::string qx_fn, std::string qy_fn, std::string mask_fn, std::string back_fn );
						
	void setBackgroundWeight( double weight );
	double backgroundWeight();
	void setOutputDirectory( std::string outdir );
	std::string outputDirectory();
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
	double p_back_weight;
	std::string p_out_dir;
	int p_alg;
	bool p_flag_subtract_background;				// subtract background?
	bool p_flag_single_correlation_output;		// write every single correlation?
	bool p_flag_use_mask;
};

#endif
