
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <string>
using std::string;

#include <sstream>		//strings
#include <fstream>		//files
#include <vector>

#include "crosscorrelator.h"
#include "arrayclasses.h"
#include "arraydataIO.h"

#include "analyze.h"

void usage(){
	cout << "====================usage============================" << endl;
	cout << " -l<filelist>   : input from a file list " << endl;
	cout << " -o<directory>  : specify output directory" << endl;
	cout << " -s             : output of separate correlation for each individual image" << endl;
	cout << " -b<filename>   : subtract background file " << endl;
	cout << " -m<filename>   : use mask file " << endl;
	cout << " -B<factor>     : weighting of background" << endl;
	cout << " -a<algorithm>  : (1) direct calculation, (2) FFT approach" << endl;
	cout << "=====================================================" << endl;
}

string argToString( char * const argument ){
	//fill argument after the -<option letter> into a string
	std::ostringstream osst;
	int j = 2;
	while ( argument[j] != '\0'){
		osst << argument[j];
		j++;
	}
	string retstring = osst.str();
	return retstring;
}

int main (int argc, char * const argv[]) {

    cout << "Hello, World!\n";

//	cout << "argc = " << argc << endl;
//	for (int i = 0; i < argc; i++){
//		cout << "argv[" << i << "] = '" << argv[i] << "'" << endl;
//	}

	Analyzer *ana = new Analyzer;
	std::vector<string> files;								// vector to be filled with individual files to process
				
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
			
				if (argv[i][1] == 'l'){
				string list_fn = argToString(argv[i]);
				
				//read from file list
				std::ifstream fin;
				fin.open(list_fn.c_str());					
				if( fin.is_open() ){
					string line;
					while (fin.good()) {
						getline(fin, line);
						if (line != ""){
							files.push_back(line);
						}
					}
				}else{
					cerr << "Could not open file '" << list_fn << "'" << endl;
					exit( 1 );
				}
				fin.close();
				

				cout << "--> using file list in " << list_fn << endl;

			}else if (argv[i][1] == 'a'){
				string a = argToString(argv[i]);
				if (a == "1"){
					ana->setAlg( 1 );
				}else if (a == "2"){
					ana->setAlg( 2 );
				}else if (a == "3"){
					ana->setAlg( 3 );
				}else if (a == "4"){
					ana->setAlg( 4 );
				}else{
					cout << "Algorithm unknown." << endl;
					exit(2);
				}
				cout << "--> using algorithm " << a << endl;	
				
			}else if (argv[i][1] == 'o'){
				string outdir = argToString(argv[i]);
				ana->setOutputDirectory( outdir );
				cout << "--> using output directory " << outdir << endl;	
				
			}else if (argv[i][1] == 's'){
				ana->flag_single_correlation_output = true;
				cout << "--> using single image output " << endl;
				
			} else if (argv[i][1] == 'b') {				
				string back_fn = argToString(argv[i]);
				
				array2D *back = new array2D;
				arraydataIO *io = new arraydataIO;
				io->readFromEDF( back_fn, back);
				ana->setBackground( back );
				ana->flag_subtract_background = true;				
				delete io;
				delete back;
				cout << "--> using background file " << back_fn << endl;
				
			} else if (argv[i][1] == 'B') {
				double weight = 1;
				std::stringstream sst;
				int j = 2;
				while ( argv[i][j] != '\0' ){				
					sst << argv[i][j];
					j++;
				}
				sst >> weight;
				ana->setBackgroundWeight( weight );
				cout << "--> using background weighting factor " << weight << endl;	
			
			}else if (argv[i][1] == 'm') {				
				string mask_fn = argToString(argv[i]);
				
				array2D *mask = new array2D;
				arraydataIO *io = new arraydataIO;
				io->readFromEDF( mask_fn, mask);
				ana->setMask( mask );			
				delete io;
				delete mask;
				cout << "--> using mask file " << mask_fn << endl;
				
			} else {
				cout << "-" << argv[i][1] << " is not a valid option." << endl;
				usage();
				return 2;
			}
		}//end for i
	}//end if
	
	
	//define set of variables to pass to the processFiles function, should probabaly go into an ini file or so at some point
	double shiftX = -7;
	double shiftY = 1;
	int num_phi = 2048;
	int num_q = 250;
	double start_q = 0;
	double stop_q = num_q;
	int LUTx = 487;				//pilatus detector x and y values
	int LUTy = 619;

	ana->processFiles( files, shiftX, shiftY, num_phi, num_q, start_q, stop_q, LUTx, LUTy);

    delete ana;
	
    cout << "Goodbye, World" << endl;
    return 0;
}



