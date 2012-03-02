//
//  tester.cpp
//  TestSuite
//
//  Created by Feldkamp on 8/23/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "tester.h"

#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <sstream>
#include <cmath>
#include <ctime>

#include <boost/array.hpp> 
using boost::array;

#include "crosscorrelator.h"
#include "arraydataIO.h"
#include "fouriertransformer.h"


//=================================================================================================
Tester::Tester(){
	const char * homepath = getenv("HOME");
	cout << "Initializing Test environment with path '" << homepath << "'" << endl;
	setBase(homepath);
}

Tester::Tester( string base ){
	setBase(base);
}

Tester::~Tester(){
}

void Tester::setBase( string base ){
	const char lastchar = base.at( base.size()-1 );
	if( lastchar != '/' ){		//if last character is not a slash, append one
		base += '/';
	}
	p_base = base;
}

string Tester::base(){
	return p_base;
}

//-------------------------------------------------------- -t1
int Tester::testCrossCorrelator( int alg, int testpattern ){						
    cout << "testCrossCorrelator(" << alg << ")" << endl;
 
	
	arraydataIO *io = new arraydataIO;   

	int sizex = 500;
	int sizey = 500;
	int nphi = 256; 
	int nq = 200;
	
	array2D<double> *dataArray = new array2D<double>(sizex, sizey);
	array2D<double> *qxArray = new array2D<double>(sizex, sizey);
	array2D<double> *qyArray = new array2D<double>(sizex, sizey);
	io->generateTestPattern(dataArray, testpattern);
	qxArray->range(-10, 10);
	qyArray->range(-10, 10);
	
    CrossCorrelator *cc = new CrossCorrelator( dataArray, qxArray, qyArray, nphi, nq );

    cc->setOutputdir( base() );
    cc->setDebug(1);							//TURN ON DEBUGGING FOR NOW --> a lot of output
    
    switch (alg) {
        case 0:
            cout << "XCCA regular" << endl;
            cc->calculatePolarCoordinates();
            cc->calculateSAXS();
            cc->calculateXCCA();	
            
        break;
        case 1:
            cout << "XCCA FAST" << endl;
			
            cc->createLookupTable(500, 500);
            double start_q = 5*cc->deltaq();
            double stop_q = cc->userQMax();

			cc->calculatePolarCoordinates_FAST(start_q, stop_q);
			cc->calculateXCCA_FAST();

			io->writeToTiff( cc->outputdir()+"polar.tif", cc->polar(), 1 );            //dump scaled output from correlation			
			io->writeToTiff( cc->outputdir()+"corr.tif", cc->autoCorr(), 1 );            //dump scaled output from correlation
			
        break;
    }

	delete io;   
	delete cc;
    return 0;
}


//-------------------------------------------------------- -t2
int Tester::testArrayClasses(){
	cout << "TESTING THE ARRAY CLASSES" << endl;

	cout << "----------array1D<double>------------" << endl;
	array1D<double> *my1Darray = new array1D<double>(5);
	my1Darray->set(0, 2);
	my1Darray->set(2, 5);
	//my1Darray->set(6, 5);//exceeds dimensions, should crash
	cout << my1Darray->getASCIIdata();
	cout << "element 2: " << my1Darray->get(2) << endl;
	cout << "element 3: " << my1Darray->get(3) << endl;	
    delete my1Darray;
	
	cout << "----------array2D<double>------------" << endl;    
	array2D<double> *my2Darray = new array2D<double>(5, 10);
	my2Darray->set(0, 0, 2);
	my2Darray->set(2, 8, 5.5);
	//my2Darray->set(6, 11, -2);//exceeds dimensions, should crash
	cout << my2Darray->getASCIIdata();
	cout << "element (2,8): " <<  my2Darray->get(2, 8) << endl;
    cout << "filling my2Darray" << endl;	
	for (int j = 0; j < my2Darray->dim2(); j++) {
		for (int i = 0; i < my2Darray->dim1(); i++) {
			my2Darray->set(i, j, j*my2Darray->dim1() + i);
		}
	}	
	cout << my2Darray->getASCIIdata();
	
	my2Darray->transpose();
	cout << "after transposing --- " << my2Darray->getASCIIdata();
	
	my2Darray->addValue(-1.5);
	cout << "after adding -1.5 --- " << my2Darray->getASCIIdata();
	my2Darray->multiplyByValue(1/2.);
	cout << "after multiplying by 1/2 --- " << my2Darray->getASCIIdata();

	array1D<double> *row = new array1D<double>();
	my2Darray->getRow(3, row);
	cout << "extracted row 3 --- " << row->getASCIIdata();
	my2Darray->setRow(4, row, 2);
	cout << "set to row 4 at index 2--- " << my2Darray->getASCIIdata();
	delete row;

	my2Darray->gradientAlongDim1(-12, 8);
	cout << "after rangeDim1(-12, 8) --- " << my2Darray->getASCIIdata();
	my2Darray->gradientAlongDim2(-12, 8);
	cout << "after rangeDim2(-12, 8) --- " << my2Darray->getASCIIdata();
	
	cout << my2Darray->getHistogramASCII( 20 );
	cout << my2Darray->getHistogramInBoundariesASCII( 20, -5, 5 );
	
	array2D<double> *mask = new array2D<double>( my2Darray->dim1(), my2Darray->dim2() );
	mask->set(1,1,1);
	mask->set(4,4,1);
	my2Darray->multiplyByArrayElementwise( mask );
	cout << "after masking: (1,1)=1 and (4,4)=1, else=0 --- " << my2Darray->getASCIIdata();
	
	delete mask;
	delete my2Darray;
	
	
    
    
//	cout << "----------array3D<double>------------" << endl;
//	array3D<double> *my3Darray = new array3D<double>(20, 10, 3);
//	my3Darray->set(0, 0, 0, 2);
//    my3Darray->set(1, 1, 1, -2);
//	my3Darray->set(2, 8, 2, 5.5);
//	my3Darray->set(6, 11, 33, -2);
//	cout << my3Darray->getASCIIdata();
//	cout << my3Darray->get(2, 8, 2) << endl;
//	cout << my3Darray->get(0, 1, 1) << endl;
//	delete my3Darray;
    
    
//    cout << "----------forging array1D<double> ---> array2D<double>----------" << endl;
//	array1D<double> *one = new array1D<double>(25);
//	one->set(0, 2);
//	one->set(2, 2);
//	one->set(6, 6);
//    one->set(24, 24);
//    array2D<double> *two = new array2D<double>(one, 5, 5);
//    cout << "one: " << one->getASCIIdata() << endl;
//    cout << "two: " << two->getASCIIdata() << endl;
//    cout << "-->deleting one" << endl;
//    delete one;
////    cout << "one: " << one->getASCIIdata() << endl;            //THIS IS SUPPOSED TO CRASH!!!
//    cout << "two: " << two->getASCIIdata() << endl;
//    cout << "-->deleting two" << endl;
//    delete two;                                                
// //   cout << "one: " << one->getASCIIdata() << endl;            //THIS IS SUPPOSED TO CRASH!!!
// //   cout << "two: " << two->getASCIIdata() << endl;            //THIS IS SUPPOSED TO CRASH!!!
 
 
 
//    cout << "----------forging array2D<double> ---> array1D<double>----------" << endl;
//    two = new array2D<double>(5, 5);
//    two->set(0,0,1);
//    two->set(1,1,7);
//    two->set(2,0,-2);
//    two->set(4,4,16);
//    one = new array1D<double>(two);
//    cout << "two: " << two->getASCIIdata() << endl;  
//    cout << "one: " << one->getASCIIdata() << endl; 
//    delete two;
////    cout << "two: " << two->getASCIIdata() << endl; 			//THIS IS SUPPOSED TO CRASH!!!
////    cout << "one: " << one->getASCIIdata() << endl; 			
//    delete one;  

    return 0;
}



//-------------------------------------------------------- -t3
int Tester::testFourierTrafo(){

    int size = 50;
    array1D<double> *f = new array1D<double>(size);
    array1D<double> *g = new array1D<double>(size);
    array1D<double> *model = new array1D<double>(size);
	
	
    //create typical 2D model
	int delta = 5;
	std::vector<int> seeds;
	seeds.push_back(6);
	seeds.push_back(24);
	seeds.push_back(32);
	seeds.push_back(3);
	seeds.push_back(22);
	seeds.push_back(1);
	for (int i = 0; i < seeds.size(); i++){
	    model->set(seeds.at(i), 1);
    	model->set(seeds.at(i)+delta, 1);
	}
	//f->addValue(1);	//put it on a pedestal
    
	//another model
	for (int i = 0; i < size; i++){
		double amp = 10.;
		double steps = 20.;	//degrees per 'i'
	    model->set(i , amp*sin( i*steps * (M_PI/180) ) );
	}
	
	
	
    cout << "--- BEFORE ALL ---" << endl;
    cout << "model: " << model->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl;
    //f->writeToASCII(base()+"/f_model.txt");
    
	double sum_model = f->calcSum();	
	double avg_model = f->calcAvg();
	
	//subtract avg SAXS intensity
	//	f->subtractValue(avg_model);
	
	
	FourierTransformer *trafo = new FourierTransformer();

    cout << "--- TRANSFORM FORWARD ---" << endl;
	f->copy( *model ); g->zeros();			//refresh arrays
	trafo->transformForward( f, g );
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl; 

    cout << "--- TRANSFORM INVERSE ---" << endl;
	f->copy( *model ); g->zeros();			//refresh arrays
	trafo->transformInverse( f, g );
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl; 
	
    cout << "--- MAGNITUDE SQUARED ---" << endl;
	f->copy( *model ); g->zeros();			//refresh arrays
	trafo->magnitudeSquared( f, g );
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl;       
    //f->writeToASCII(base()+"/f_power.txt");
    
    cout << "--- CORRELATION ---" << endl;
	f->copy( *model ); g->zeros();			//refresh arrays
    trafo->autocorrelation(f, g);
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl;      
    //f->writeToASCII(base()+"/f_corr.txt");
	
	
	//normalize according to Altarelli
	//f->multiplyByValue( 1/avg_model/avg_model);
	
  	double sum_corr = f->calcSum();
	double avg_corr = f->calcAvg();
	cout << "sum_model=" << sum_model << ", avg_model=" << avg_model << endl;
	cout << "sum_corr=" << sum_corr << ", avg_corr=" << avg_corr << endl;  
	
	delete trafo;
    delete f;
    delete g;
    
	return 0;
}


//-------------------------------------------------------- -t4
int Tester::testIO(int mode){
	cout << "TESTING THE ARRAYDATA INPUT/OUTPUT" << endl;
	
	string outfn = "test";
	string infn = "test";
	
	arraydataIO *io = new arraydataIO(2);

	//-------------------------------------------------------------writing
	cout << "--- test case 2D ---" << endl;
	array2D<double> *wmat = new array2D<double>(10, 6);
	io->generateTestPattern( wmat, 1 );		//cases 0 - 4 available
//	wmat->gradientAlongDim1(-10, 10);

	cout << wmat->getASCIIdata();
	
	
	if (mode == 0 || mode == 1){
		cout << "writing to EDF" << endl;
		io->writeToEDF(base()+outfn+".edf", wmat);
	}
	
	if (mode == 0 || mode == 2){
		cout << "writing to HDF5 (doubles)" << endl;
		io->writeToHDF5(base()+outfn+"_d.h5", wmat, 0);
		cout << "writing to HDF5 (float)" << endl;
		io->writeToHDF5(base()+outfn+"_f.h5", wmat, 1);
		cout << "writing to HDF5 (int)" << endl;
		io->writeToHDF5(base()+outfn+"_i.h5", wmat, 2);
		cout << "writing to HDF5 (int16_t)" << endl;
		io->writeToHDF5(base()+outfn+"_i16.h5", wmat, 3);
		cout << "writing to HDF5 (long)" << endl;
		io->writeToHDF5(base()+outfn+"_long.h5", wmat, 4);
	}
	
	if (mode == 0 || mode == 3){
		cout << "writing to 16-bit TIFF scaled" << endl;
		io->writeToTiff(base()+outfn+"_16scaled.tif", wmat, 1);	
		cout << "writing to 16-bit TIFF unscaled" << endl;
		io->writeToTiff(base()+outfn+"_16unscaled.tif", wmat, 0);
		cout << "writing to 32-bit TIFF" << endl;
		io->writeToTiff(base()+outfn+"_32.tif", wmat, 0, 32);
	}
	
	if (mode == 0 || mode == 4){
		cout << "writing to ASCII (1D)" << endl;
		io->writeToASCII(base()+outfn+"_1D.txt", wmat, 1);	
		cout << "writing to ASCII (2D)" << endl;
		io->writeToASCII(base()+outfn+"_2D.txt", wmat);
	}
	
	//-------------------------------------------------------------reading
	array2D<double> *rmat = new array2D<double>;

	if (mode == 0 || mode == 1){
		cout << "\nreading from EDF" << endl;
		io->readFromEDF(base()+infn+".edf", rmat);
		cout << rmat->getASCIIdata();
	}
	
	if (mode == 0 || mode == 2){
		cout << "\nreading from HDF5(double)" << endl;
		io->readFromHDF5(base()+infn+"_d.h5", rmat);
		cout << rmat->getASCIIdata();
		cout << "\nreading from HDF5(float)" << endl;
		io->readFromHDF5(base()+infn+"_f.h5", rmat);
		cout << rmat->getASCIIdata();
		cout << "\nreading from HDF5(int)" << endl;
		io->readFromHDF5(base()+infn+"_i.h5", rmat);
		cout << rmat->getASCIIdata();
		cout << "\nreading from HDF5(int16_t)" << endl;
		io->readFromHDF5(base()+infn+"_i16.h5", rmat);
		cout << rmat->getASCIIdata();
		cout << "\nreading from HDF5(long)" << endl;
		io->readFromHDF5(base()+infn+"_long.h5", rmat);
		cout << rmat->getASCIIdata();
	}
		
	if (mode == 0 || mode == 3){
		cout << "\nreading from 16-bit TIFF scaled" << endl;
		io->readFromTiff(base()+infn+"_16scaled.tif", rmat);
		cout << rmat->getASCIIdata();
		cout << "\nreading from 16-bit TIFF unscaled" << endl;
		io->readFromTiff(base()+infn+"_16unscaled.tif", rmat);
		cout << rmat->getASCIIdata();
		cout << "\nreading from 32-bit TIFF" << endl;
		io->readFromTiff(base()+infn+"_32.tif", rmat);
		cout << rmat->getASCIIdata();
	}
	
	if (mode == 0 || mode == 4){
		cout << "\nreading from ASCII (2D)" << endl;
		io->readFromASCII(base()+infn+"_2D.txt", rmat);
		cout << rmat->getASCIIdata();
	}
	
	//special testing mode.....
	if (mode == 5){
		cout << "\nreading from large file ASCII (2D)" << endl;
		io->readFromASCII("/Users/feldkamp/Work/SLAC/2011_06_Water_L357/results/kitty_out/PEDESTALSNEW/cspad-pedestals_r0006.dat", rmat);	
	}
	
	
	//-------------------------------------------------------------writing
	cout << "--- test case 1D ---" << endl;
	array1D<double> *warr = new array1D<double>(10);
	warr->range(-7, 12);
	cout << warr->getASCIIdata();

	if (mode == 10 || mode == 11) {
		cout << "writing to EDF (1D)" << endl;
		io->writeToEDF(base()+outfn+"_1D.edf", warr);
	}
	
	if (mode == 10 || mode == 12) {
		cout << "writing to HDF5 (1D)" << endl;
		io->writeToHDF5(base()+outfn+"_1D.h5", warr);	
	}

	if (mode == 10 || mode == 14) {
		cout << "writing to ASCII (1D)" << endl;	
		io->writeToASCII(base()+outfn+"_1D.txt", warr);
	}
	
	//-------------------------------------------------------------reading
	array1D<double>* rarr = new array1D<double>;
	
	if (mode == 10 || mode == 11) {
		io->readFromEDF(base()+infn+"_1D.edf", rarr);
		cout << rarr->getASCIIdata();
	}

	if (mode == 10 || mode == 12) {
		cout << "\nreading from HDF5 (1D)" << endl;
		io->readFromHDF5(base()+infn+"_1D.h5", rarr);
		cout << rarr->getASCIIdata();
	}

	if (mode == 10 || mode == 14) {
		cout << "\nreading from ASCII (1D)" << endl;
		io->readFromASCII(base()+infn+"_1D.txt", rarr);
		cout << rarr->getASCIIdata();
	}

	delete wmat;
	delete warr;
	delete rmat;	
	delete rarr;
	delete io;
	
	return 0;
}


//-------------------------------------------------------- -t5
int Tester::testDataTypes(){
	cout << "--- TESTING DATA TYPES ON THIS SYSTEM ---" << endl;
	cout << "sizeof(int) =                " << sizeof(int) << endl;
	cout << "sizeof(short int) =          " << sizeof(short int) << endl;
	cout << "sizeof(long int) =           " << sizeof(long int) << endl;
	cout << "sizeof(unsigned int) =       " << sizeof(unsigned int) << endl;	
	cout << "sizeof(unsigned short int) = " << sizeof(unsigned short int) << endl;	
	cout << "sizeof(unsigned long int)  = " << sizeof(unsigned long int) << endl;		
	cout << "sizeof(int16_t)  =           " << sizeof(int16_t) << endl;
	cout << "sizeof(int32_t)  =           " << sizeof(int32_t) << endl;
	cout << "sizeof(float) =              " << sizeof(float) << endl;
	cout << "sizeof(double) =             " << sizeof(double) << endl;
	cout << "sizeof(bool) =               " << sizeof(bool) << endl;
	cout << "-----------------------------------------" << endl;
	return 0;
}




//-------------------------------------------------------- -t6
int Tester::testArraySpeed(){ 
	int size1D = 1000*1000;
	unsigned int dim1 = 10000;
	unsigned int dim2 = 10000;
	
	clock_t c_start, c_bst, c_rry, c_vctr, c_2D_a, c_2D_b;
	
	c_start = clock();
	subtestBoost();
	c_bst = clock() - c_start;

	c_start = clock();
	subtestArraydata(size1D);
	c_rry = clock() - c_start;

	c_start = clock();	
	subtestVector(size1D);
	c_vctr = clock() - c_start;

	//specific 2D tests
	c_start = clock();	
	subtest2D_a(dim1,dim2);
	c_2D_a = clock() - c_start;

	c_start = clock();	
	subtest2D_b(dim1,dim2);
	c_2D_b = clock() - c_start;
	
	cout << "1D boost  " << c_bst << ", " << c_bst/double(CLOCKS_PER_SEC)*1000 << " ms" << endl;
	cout << "1D array  " << c_rry << ", " << c_rry/double(CLOCKS_PER_SEC)*1000 << " ms" << endl;
	cout << "1D vector " << c_vctr << ", " << c_vctr/double(CLOCKS_PER_SEC)*1000 << " ms" << endl;
	cout << "2D_a test " << c_2D_a << ", " << c_2D_a/double(CLOCKS_PER_SEC)*1000 << " ms" << endl;
	cout << "2D_b test " << c_2D_b << ", " << c_2D_b/double(CLOCKS_PER_SEC)*1000 << " ms" << endl;
	cout << "ratio b/a: " << c_2D_a/(double)c_2D_b << endl;
	
	return 0;
} 


void Tester::subtestBoost(){
	cout << "subtestBoost" << endl;
	const int size = 1000;
	array<double, size> a;
	for (unsigned int i = 0; i < a.size(); i++){
		a.at(i) = 0;
	}
}

void Tester::subtestArraydata(unsigned int size){
	cout << "subtestArraydata" << endl;
	arraydata<double> a = arraydata<double>(size);
	for (unsigned int i = 0; i < a.size(); i++){
		a.set_atIndex(i, double(i));
		double value = a.get_atIndex(i);
//		if (!(i%1000)){
//			cout << value << " " << std::flush;
//			if (!(i%1000000))
//				cout << endl;
//		}
	}
}

void Tester::subtestVector(unsigned int size){
	cout << "subtestVector" << endl;
	vector<double> v(size, 0);

	for (unsigned int i = 0; i < v.size(); i++){
		v.at(i) = double(i);
		double value = v.at(i);
//		if (!(i%1000)){
//			cout << value << " " << std::flush;
//			if (!(i%1000000))
//				cout << endl;
//		}
	}
}

void Tester::subtest2D_a(unsigned int dim1, unsigned int dim2){
	cout << "subtest2D_a: i slow, j fast variable in loop" << endl;
	array2D<double> img(dim1, dim2);
	for (int i = 0; i < img.dim1(); i++){			// i
		for (int j = 0; j < img.dim2(); j++){		// j
			img.set(i, j, 1);
		}
	}
}

void Tester::subtest2D_b(unsigned int dim1, unsigned int dim2){
	cout << "subtest2D_b: j slow, i fast variable in loop" << endl;
	array2D<double> img(dim1, dim2);
	for (int j = 0; j < img.dim2(); j++){			// j
		for (int i = 0; i < img.dim1(); i++){		// i
			img.set(i, j, 1);
		}
	}
}



