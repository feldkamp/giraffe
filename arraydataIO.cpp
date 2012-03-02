//
//  arraydataIO.cpp
//  xcca_commandline
//
//  Created by Feldkamp on 7/3/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "arraydataIO.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <sstream>				//string steaming

#include <fstream>				//file streaming

#include <iterator>

#include <cmath>

//-------------------------------------------------------------getExt

//-------------------------------------------------------------
std::string getExt(std::string filename){
	char separator = '.';
	//find separator dot from the end of the string
	size_t pos = filename.rfind( separator );
	//get extension without the '.'
	string ext = filename.substr( pos+1, filename.length() );
	return ext;
}


arraydataIO::arraydataIO( int verbose )
	: p_transposeForIO(true)	//should be 'true', generally. see note in header
	, p_verbose(verbose)
{
}



arraydataIO::~arraydataIO(){

}

int arraydataIO::readFromFile( std::string filename, array1D<double> *&dest) const{
	string ext = getExt( filename );
	int fail = 0;
	
	// read 2D 'raw' image
	if (ext == "edf" || ext == "EDF" || ext == "bin"){			
		fail = readFromEDF( filename, dest );
	}else if (ext == "h5" || ext == "hdf5" || ext == "H5" || ext == "HDF5"){		
		fail = readFromHDF5( filename, dest );
	}else if (ext == "txt" || ext == "TXT" || ext == "dat"){		
		fail = readFromASCII( filename, dest );
	}else{
		cerr << "Error in arraydataIO::readFromFile(1D)! Extension '" << ext << "' found in '" << filename << "' is not valid. " << endl;
		cerr << "valid options include 'edf', 'h5' and 'txt/dat'" << endl;
		fail = 100;
	}
	return fail;
}

int arraydataIO::readFromFile( std::string filename, array2D<double> *&dest) const{
	string ext = getExt( filename );
	int fail = 0;
	
	// read 2D 'raw' image
	if (ext == "edf" || ext == "EDF" || ext == "bin"){			
		fail = readFromEDF( filename, dest );
	}else if (ext == "h5" || ext == "hdf5" || ext == "H5" || ext == "HDF5"){		
		fail = readFromHDF5( filename, dest );
	}else if (ext == "tif" || ext == "tiff" || ext == "TIF" || ext == "TIFF"){		
		fail = readFromTiff( filename, dest );
	}else if (ext == "txt" || ext == "TXT" || ext == "dat"){		
		fail = readFromASCII( filename, dest );
	}else{
		cerr << "Error in arraydataIO::readFromFile(2D)! Extension '" << ext << "' found in '" << filename << "' is not valid. " << endl;
		cerr << "valid options include 'edf', 'h5', 'tif', and 'txt/dat'" << endl;
		fail = 100;
	}
	return fail;
}

int arraydataIO::writeToFile( std::string filename, array1D<double> *src) const{
	string ext = getExt( filename );
	int fail = 0;
	
	if (ext == "edf" || ext == "EDF" || ext == "bin"){			
		fail = writeToEDF( filename, src );
	}else if (ext == "h5" || ext == "hdf5" || ext == "H5" || ext == "HDF5"){		
		fail = writeToHDF5( filename, src );
	}else if (ext == "txt" || ext == "TXT" || ext == "dat"){		
		fail = writeToASCII( filename, src );
	}else{
		cerr << "Error in arraydataIO::writeToFile(1D)! Extension '" << ext << "' found in '" << filename << "' is not valid. " << endl;
		cerr << "valid options include 'edf', 'h5' and 'txt/dat'" << endl;
		fail = 100;
	}
	return fail;
}

int arraydataIO::writeToFile( std::string filename, array2D<double> *src) const{
	string ext = getExt( filename );
	int fail = 0;

	if (ext == "edf" || ext == "EDF" || ext == "bin"){			
		fail = writeToEDF( filename, src );
	}else if (ext == "h5" || ext == "hdf5" || ext == "H5" || ext == "HDF5"){		
		fail = writeToHDF5( filename, src );
	}else if (ext == "tif" || ext == "tiff" || ext == "TIF" || ext == "TIFF"){
		fail = writeToTiff( filename, src );
	}else if (ext == "txt" || ext == "TXT" || ext == "dat"){		
		fail = writeToASCII( filename, src );
	}else{
		cerr << "Error in arraydataIO::writeToFile(2D)! Extension '" << ext << "' found in '" << filename << "' is not valid. " << endl;
		cerr << "valid options include 'edf', 'h5', 'tif', and 'txt/dat'" << endl;
		fail = 100;
	}
	return fail;
}



//****************************************************************************************************************
//****************************************************************************************************************
//*** EDF                                                                                                      ***
//****************************************************************************************************************
//****************************************************************************************************************
#ifdef ARRAYDATAIO_EDF
	#include "edf.h"				// needed to read/write EDF files (ESRF Data Format)


	//-------------------------------------------------------------- readFromEDF (1D)
	// implementation borrowed following vectordata/matrixdata::in() 
	// of the tomo/matlib package (see http://xray-lens.de )
	//--------------------------------------------------------------
	int arraydataIO::readFromEDF( string filename, array1D<double> *&dest ) const{
		ns_edf::edf *file = new ns_edf::edf;
		if (file->read_header(filename)){
			cerr << "Error in arraydataIO::readFromEDF! Could not read EDF file '" << filename << "'." << endl;
			return 1;
		}
		
		if ( file->get_Dim_y() > 1 )
			cerr << "Warning in arraydataIO::readFromEDF! EDF-file is not 1D!" << endl;			
		int dim1 = (int)file->get_Dim_x();
		
		if (verbose()){
			cout << "Reading '" << filename << "' (EDF file type " << file->get_FileType_str() << ")"
			<< ", dimension " << dim1 << "" << endl;
		}
		double *temp = new double[dim1];
		
		int fail = file->read_data( temp, filename );
		if ( fail ){
			cerr << "ERROR. In arraydataIO::readFromEDF - Could not load EDF data file '" << filename << "'." << endl;
			return 2;
		}
		
		//feed back into dest
		delete dest;
		dest = new array1D<double>( temp, dim1 );
		
		delete []temp;
		delete file;
		return 0; 
	}
		
	//-------------------------------------------------------------- readFromEDF (2D)
	//
	//
	//--------------------------------------------------------------------------
	int arraydataIO::readFromEDF( string filename, array2D<double> *&dest ) const{
		ns_edf::edf *file = new ns_edf::edf;

		if (file->read_header(filename)){
			cerr << "Error in arraydataIO::readFromEDF! Could not read EDF file '" << filename << "'." << endl;
			return 1;
		}
		
		int dim1 = (int)file->get_Dim_x();
		int dim2 = (int)file->get_Dim_y();
		
		if(verbose()){
			cout << "Reading '" << filename << "' (EDF file type " << file->get_FileType_str() << ")"
			<< ", dimensions (" << dim1 << ", " << dim2 << ")" << endl;
		}	
		double *temp = new double[dim1*dim2];
		
		int fail = file->read_data( temp, filename );
		if ( fail ){
			cerr << "ERROR. In arraydataIO::readFromEDF - Could not load EDF data file '" << filename << "'." << endl;
			return 2;
		}
		
		//feed back into dest
		delete dest;
		dest = new array2D<double>( temp, dim1, dim2 );
		//dest->arraydata::copy( temp, dim1*dim2);
		
		delete []temp;
		
		if (transposeForIO()){
			dest->transpose();
		}

		delete file;
		return 0; 
	}

	//-------------------------------------------------------------- writeToEDF (generic)
	// non-class function
	//--------------------------------------------------------------
	int writeToEDF_generic( ns_edf::edf *file, arraydata<double> *src ) {
	
		bool flipByteOrder = false;
		
		// write scaled data?
		// SF_UNSCALED = 0:  binary values written represent actual data values (a range from 0 to 65535 is possible)
		// SF_SCALED	= 1: scaling max and min are written to header, data values always span full range from 0 to 65535 (2^16-1)
		// SF_RETAIN	= 2: if the scaling factor was read from an input file, keep it what it was
		ns_edf::scaled_t scaleOut = ns_edf::SF_SCALED;
		
		file->set_ScalingMin( src->calcMin() );
		file->set_ScalingMax( src->calcMax() );
		file->set_ScaledFlag( scaleOut );
		
		double *buf = new double[src->size()];
		for (int i = 0; i < src->size(); i++){
			buf[i] = src->get_atIndex(i);
		}
		return file->write( buf, flipByteOrder );	
	}

	//-------------------------------------------------------------- writeToEDF (1D)
	//
	//--------------------------------------------------------------
	int arraydataIO::writeToEDF( string filename, array1D<double> *src ) const {
		if ( !src ){
			cerr << "ERROR. In arraydataIO::writeToEDF(1D). No source data. Could not write to file " << filename << endl;
			return 1;
		}

		ns_edf::edf *file = new ns_edf::edf;
		
		file->set_Dim_x( src->dim1() );
		file->set_Dim_y( 1 );
		file->set_FilenameOut( filename );
		int retval = writeToEDF_generic( file, src );
		
		delete file;
		return retval;
	}

	//-------------------------------------------------------------- writeToEDF (2D)
	//
	//--------------------------------------------------------------
	int arraydataIO::writeToEDF( string filename, array2D<double> *src ) const {
		if ( !src ){
			cerr << "ERROR. In arraydataIO::writeToEDF(2D). No source data. Could not write to file " << filename << endl;
			return 1;
		}

		if (transposeForIO()){//transpose before writing
			src->transpose();
		}
		ns_edf::edf *file = new ns_edf::edf;
		
		file->set_Dim_x( src->dim1() );
		file->set_Dim_y( src->dim2() );
		file->set_FilenameOut( filename );		
		int retval = writeToEDF_generic( file, src );
		
		delete file;
		
		if (transposeForIO()){//transpose back before returning
			src->transpose();
		}
		return retval;
	}

#else
	int arraydataIO::readFromEDF( string filename, array1D<double> *&dest ) const{
		cerr << "======== WARNING! Dummy function. Nothing read from EDF. ======== " << endl; 
		return 1;
	}
	int arraydataIO::readFromEDF( string filename, array2D<double> *&dest ) const{
		cerr << "======== WARNING! Dummy function. Nothing read from EDF. ======== " << endl; 
		return 1;
	}
	int arraydataIO::writeToEDF( string filename, array1D<double> *src ) const {
		cerr << "======== WARNING! Dummy function. Nothing written to EDF.======== " << endl; 
		return 1;
	}
	int arraydataIO::writeToEDF( string filename, array2D<double> *src ) const {
		cerr << "======== WARNING! Dummy function. Nothing written to EDF.======== " << endl; 
		return 1;
	}
#endif









//****************************************************************************************************************
//****************************************************************************************************************
//*** TIFF                                                                                                     ***
//****************************************************************************************************************
//****************************************************************************************************************
#ifdef ARRAYDATAIO_TIFF
	#include <tiffio.h>				// needed to read/write tiff files (Tagged Image File Format)
	
	//-------------------------------------------------------------- readFromTiff
	// (needs TIFFLIB installation)
	//
	// implementation borrowed following matrixdata::tiffin 
	// of the tomo/matlib package (see http://xray-lens.de )
	//--------------------------------------------------------------------------
	int arraydataIO::readFromTiff( string filename, array2D<double> *&dest ) const {
		int retval = 0;
		TIFF *tiff= TIFFOpen( filename.c_str(), "r" );

		if(tiff){
			uint32 width, imagelength;
			uint32 datalength;
			uint32 datatype = 0;
			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &imagelength);
			TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &datalength);
			TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &datatype);
							
			//allocate data to write image to		
			delete dest;
			dest = new array2D<double>(width, imagelength);

			if(verbose()){
				cout << "Reading '" << filename << "' (TIFF file)"
				<< ", dimensions (" << dest->dim1() << ", " << dest->dim2() << "), " << datalength << " bit." << endl;
			}
						
			switch ( datalength) {
			case 32:

				if( datatype == SAMPLEFORMAT_UINT )
				{
					uint32* ibuf;
					ibuf = (uint32*) _TIFFmalloc(TIFFScanlineSize(tiff));				

					for (uint32 row = 0; row < imagelength; row++)
					{
						TIFFReadScanline(tiff, ibuf, row);
						for (uint16 i = 0; i < width; i++)
							dest->set(i, row, ibuf[i]);
					}
					_TIFFfree(ibuf);
				}
				else if( datatype == SAMPLEFORMAT_IEEEFP )
				{
					float* fbuf;
					fbuf = (float*) _TIFFmalloc(TIFFScanlineSize(tiff));
					for (uint32 row = 0; row < imagelength; row++)
					{
						TIFFReadScanline(tiff, fbuf, row);
						for (uint16 i = 0; i < width; i++)
							dest->set(i, row, fbuf[i]);
					}
					_TIFFfree(fbuf);
				}else{
					cerr << "arraydataIO::readFromTiff. Data format " << datatype << " not implemented." << endl;
				}

				break;
			case 16:
				uint16* buf;

				buf = (uint16*) _TIFFmalloc(TIFFScanlineSize(tiff));

				for (uint32 row = 0; row < imagelength; row++)
				{
					TIFFReadScanline(tiff, buf, row);
					for (uint16 i = 0; i < width; i++)
					{
						dest->set(i, row, buf[i]);
					}
				}

				_TIFFfree(buf);
				break;
			default:
				int npixels;
				uint32* raster;

				npixels = width * imagelength;
				raster = (uint32*) _TIFFmalloc( npixels * sizeof(uint32) );
				if (raster != NULL){
					if ( TIFFReadRGBAImage(tiff, width, imagelength, raster, 0) ){
						for(uint32 row= 0; row < imagelength; row++){		//y
							for(uint16 i= 0; i < width; i++){				//x	
	//							dest->set( i, j, 1.* (raster[j*width+i] & 0x000000ff));			// max value 255
								dest->set( i, row, 1.* (raster[row*width+i] & 0x0000ffff));			// max value 65535, force to double
							}
						}
					}
					_TIFFfree(raster);
					dest->fliplr();
				}else{
					cerr << "Error in arraydataIO::readFromTiff. Could not read image (TIFFReadRGBAImage failed)." << endl;
				}
			}//switch	
			
			if (transposeForIO()){
				dest->transpose();
			}
			
			TIFFClose(tiff);
		} else { //tiff could not be opened
			retval = 1;
		}
		
		return retval;
	}





	//-------------------------------------------------------------- writeToTiff
	// (needs TIFFLIB installation)
	// ATTENTION: libtiff is not explicitly thread-safe
	//			be sure to use mutexes when writing to disk
	//
	// scaleFlag 0: direct output of values, 
	//				cut off values larger than 2^16-1 = 65535
	// scaleFlag 1: scaled output to full tiff scale of 0 to 65535
	//
	// implementation borrowed following matrixdata::tiff16out 
	// of the tomo/matlib package (see http://xray-lens.de )
	//--------------------------------------------------------------------------
	int arraydataIO::writeToTiff( string filename, array2D<double> *src, int scaleFlag ) const {
		if ( !src ){
			cerr << "ERROR. In arraydataIO::writeToTiff. No source data. Could not write to file " << filename << endl;
			return 1;
		}
		
		if (src->size() == 0) {
			cerr << "Error in writeToTiff! Array size is zero." << endl;
			return 2;
		}
		
		TIFF *out = TIFFOpen(filename.c_str() ,"w");
		if(out)
		{
			if (transposeForIO()){
				src->transpose();
			}
			const uint32 width = (uint32) src->dim1();
			const uint32 height = (uint32) src->dim2();
			
			const double MaxValue = src->calcMax();
			const double MinValue = src->calcMin();
			const double MaxMinRange = MaxValue - MinValue;
			
			uint16 *tifdata = new uint16[width];
			
			// define tif-parameters
			const uint16 spp = 1;      // Samples per pixel
			const uint16 bpp = 16;     // Bits per sample
			const uint16 photo = PHOTOMETRIC_MINISBLACK;
			
			// set the width of the image
			TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width/spp);  
			//TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width * bpp * 8); 
			
			// set the height of the image
			TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);  
			
			// set the size of the channels   
			TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp); 
			
			// set number of channels per pixel
			TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
			
			// set the origin of the image
			TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    
			
			TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
			TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);

			//write to file
			int too_low_count = 0;
			int too_high_count = 0;
			for (uint32 i = 0; i < height; i++)
			{
				for (unsigned int j = 0; j < width; j++)
				{
					if (scaleFlag) {										//scale image to full range
						tifdata[j] = (uint16) floor( 65535.* (src->get(j,i)-MinValue) / MaxMinRange );
					} else {										
						if( src->get(j,i) < 0 ) {
							tifdata[j] = 0;									// cut off values smaller than 0
							too_low_count++;
						}
						else if ( src->get(j,i) > 65535 ) {
							tifdata[j] = 65535;								// ...and values larger than 2^16-1
							too_high_count++;
						}
						else{ 
							tifdata[j] = (uint16) floor(src->get(j,i));		// ... force everything else to be integer
						}
					}
				}
				
				TIFFWriteScanline(out, tifdata, i, 0);
			}
			
			
			
			if (verbose()){
				cout << "Tiff image '" << filename << "' written to disc." << endl;
				if (scaleFlag) cout << "Scaled output. "; else cout << "Unscaled output. ";
				cout << "data min: " << MinValue << ", max: " << MaxValue << ", range: " << MaxMinRange << endl;
				
				if (too_low_count){
					cout << too_low_count << " values lower than zero in the original data have been cut off." << endl;
				}
				if (too_high_count){
					cout << too_high_count << " values higher than 65535 in the original data have been cut off." << endl;
				}
			}
			
			delete[] tifdata;
			TIFFClose(out);
			
			if (transposeForIO()){//transpose back before returning
				src->transpose();
			}
			return 0;
		}else{
			cerr << "Error in array2D<double>::writeToTiff! Could not open '" << filename << "' for writing!" << endl;
			return 1;	
		}
	}

#else
	int arraydataIO::readFromTiff( string filename, array2D<double> *&dest ) const {
		cerr << "======== WARNING! Dummy function. Nothing read from Tiff. ======== " << endl; 
		return 1;
	}
	int arraydataIO::writeToTiff( string filename, array2D<double> *src, int scaleFlag ) const {
		cerr << "======== WARNING! Dummy function. Nothing written to Tiff. ======== " << endl; 
		return 1;
	}
#endif






//****************************************************************************************************************
//****************************************************************************************************************
//*** HDF5                                                                                                     ***
//****************************************************************************************************************
//****************************************************************************************************************

#ifdef ARRAYDATAIO_HDF5
	#include <hdf5.h>				// needed to read/write HDF5 files (Tagged Image File Format)
	
	//-------------------------------------------------------------- readFromHDF5 (generic case)
	//
	//-------------------------------------------------------------- 
	int readFromHDF5_generic( string filename, arraydata<double> *&dest, int &rank, hsize_t *&dims, int verbose ) {
		
		if (verbose){
			cout << "Reading " << filename << " (HDF5 file)" << endl;
		}	
		
		std::ostringstream info;
		
		string dataset_name = "data/data";
		unsigned int access_mode = H5F_ACC_RDONLY;			// file mode, read-only, other option: H5F_ACC_RDWR
		hid_t access_prp = H5P_DEFAULT;						// file access property list
		
		hid_t fileID = H5Fopen(filename.c_str(), access_mode, access_prp);			// open file
		if (fileID <= 0){
			cerr << "ERROR in arraydataIO::readFromHDF5. Could not read HDF5 file '" << filename << "'!" << endl;
			return 1;
		}
		hid_t datasetID = H5Dopen(fileID, dataset_name.c_str(), access_prp);		// open dataset in that file


		hid_t memTypeID = H5T_NATIVE_INT;		//datatype in memory (not necessarily in the file)
		hid_t memSpaceID = H5S_ALL;				//dataspace in memory, H5S_ALL is default --> the whole dataspace in memory will be used for I/O
		hid_t fileSpaceID = H5S_ALL;			//dataspace in file, H5S_ALL is default --> the whole dataspace in the fill is selected for I/O
		hid_t xfer_prp = H5P_DEFAULT;			//transfer property list, H5P_DEFAULT is default
				
		hid_t typeID = H5Dget_type(datasetID);
		H5T_class_t varClass = H5Tget_class(typeID);

		hid_t dataspace = H5Dget_space(datasetID);    /* dataspace handle */
		rank = H5Sget_simple_extent_ndims(dataspace);
		delete dims;
		dims = new hsize_t[rank];
		int status_n = H5Sget_simple_extent_dims(dataspace, dims, NULL);
		
		info << "rank " << rank << ". " << flush;
		int n, rows, cols = 0;
		if (rank == 1){
			rows = (int)(dims[0]);
			n = rows;
			info << "rows = " << rows << ". " << flush;
		}else if (rank == 2){
			rows = (int)(dims[0]);
			cols = (int)(dims[1]);
			n = rows * cols;
			info << "rows,cols = " << rows << "," << cols << ". " << flush;
		}else{
			cerr << "--> can't handle that rank. Aborting." << endl;
			return 1;
		}
		info << "(" << n << "). status: " << status_n << ". " << flush;

		H5T_order_t order = H5Tget_order(typeID);
		if (order == H5T_ORDER_LE){
			info << "LE. " << flush;	//little endian byte order
		}else if (order == H5T_ORDER_BE){
			info << "BE. " << flush;	//big endian byte order
		}
		size_t varSize = H5Tget_size(typeID);

		//allocate data to hold what comes out of the file
		void *buffer = NULL;
		if (varClass == H5T_FLOAT){
			info << "H5type FLOAT, size " << (int)varSize << " " << flush;
			if (varSize == sizeof(double)){
				info << "(double) " << flush;
				buffer = new double[n];
				memTypeID = H5T_NATIVE_DOUBLE;
			}else if (varSize == sizeof(float)){
				info << "(float) " << flush;
				buffer = new float[n];
				memTypeID = H5T_NATIVE_FLOAT;
			}
		}else if (varClass == H5T_INTEGER){
			info << "H5type INTEGER, size " << (int)varSize << " " << flush;
			if (varSize == sizeof(int)){
				info << "(int) " << flush;
				buffer = new int[n];
				memTypeID = H5T_NATIVE_INT;
			}else if (varSize == sizeof(int16_t)){
				info << "(int16_t) " << flush;
				buffer = new int16_t[n];
				memTypeID = H5T_STD_I16LE;
			}else if (varSize == sizeof(long int)){
				info << "(long) " << flush;
				buffer = new long int[n];
				memTypeID = H5T_NATIVE_LONG;
			}
			
			H5T_sign_t sign = H5Tget_sign(typeID);
			if (sign == H5T_SGN_NONE){
				info << "unsigned " << flush;
			}else if (sign == H5T_SGN_2){
				info << "signed " << flush;
			}
		}else{
			cerr << "Error in arraydataIO::readFromHDF5. Unknown HDF5 class type '" << varClass << "'" << endl;
			return 1;
		}
		
		if (!buffer){
			cerr << "Error in readFromHDF5, could not allocate read buffer of size " << n << "." << endl;
			return 2;
		}

		// read the data using the default properties.
		//hid_t status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);
		hid_t status = H5Dread (datasetID, memTypeID, memSpaceID, fileSpaceID, xfer_prp, buffer);
		
		// close and release resources.
		status = H5Dclose(datasetID);
		status = H5Fclose(fileID);
				
		//copy data back to dest, then delete buffer memory
		//-----ATTENTION-------
		//for some reason, the dereferencing of anything else than 'int*' seems to fail for floating point types
		//this throws away all floating point accuracy in the data.......
		if (memTypeID == H5T_NATIVE_DOUBLE) {
			delete dest;
			dest = new arraydata<double>( (double*)buffer, n );
			delete[] ((double*)buffer);
		} else if (memTypeID == H5T_NATIVE_FLOAT) {
			delete dest;
			dest = new arraydata<double>( (float*)buffer, n );
			delete[] ((float*)buffer);
		} else if (memTypeID == H5T_NATIVE_INT) { 
			delete dest;
			dest = new arraydata<double>( (int*)buffer, n );
			delete[] ((int*)buffer);
		} else if (memTypeID == H5T_STD_I16LE) { 
			delete dest;
			dest = new arraydata<double>( (int16_t*)buffer, n );
			delete[] ((int16_t*)buffer);
		} else if (memTypeID == H5T_NATIVE_LONG) { 
			delete dest;
			dest = new arraydata<double>( (long*)buffer, n );
			delete[] ((long*)buffer);
		} else { 
			cerr << "Error in arraydataIO::readFromHDF5. Type not found. " << endl;
			return 1;
		}
		
		if (verbose>1){
			cout << info.str() << endl;
		}
		return 0;
	}

	//-------------------------------------------------------------- readFromHDF5 (return array1D<double> object)
	int arraydataIO::readFromHDF5( string filename, array1D<double> *&dest ) const {
		int img_rank = 0;
		hsize_t *dims = 0;
		arraydata<double> *readarray = 0;
		int fail = readFromHDF5_generic( filename, readarray, img_rank, dims, verbose() );
		if (!fail){
			delete dest;
			dest = new array1D<double>( readarray );
		}
		delete readarray;
		return fail;		
	}

	//-------------------------------------------------------------- readFromHDF5 (return array2D<double> object)
	int arraydataIO::readFromHDF5( string filename, array2D<double> *&dest ) const {
		int img_rank = 0;
		hsize_t *dims = 0;
		arraydata<double> *readarray = 0;
		int fail = readFromHDF5_generic( filename, readarray, img_rank, dims, verbose() );
		if (!fail){
			int cols = (int)dims[0];
			int rows = (int)dims[1];
			if (img_rank == 1){ cols = 1; }
			convertToArray2D( readarray, dest, rows, cols );
			if (verbose()>2){
				cout << "dest: " << dest << ", dims " << dest->dim1() << " x " << dest->dim2() 
					<< ", rows:" << rows << ", cols:" << cols << endl;
				cout << dest->getASCIIdata();
			}
			
			if (transposeForIO()){
				dest->transpose();
			}
		}
		delete readarray;

		return fail;
	}	

	//-------------------------------------------------------------- writeToHDF5
	// dataType = 0 --> write as doubles (H5T_NATIVE_DOUBLE, default)
	// dataType = 1 --> write as float   (H5T_NATIVE_FLOAT)
	// dataType = 2 --> write as int     (H5T_NATIVE_INT)
	// dataType = 3 --> write as int16_t (H5T_STD_I16LE)
	// dataType = 4 --> write as long    (H5T_NATIVE_LONG)
	//-------------------------------------------------------------- 
	//-------------------------------------------------------------- writeToHDF5 (generic case)
	int writeToHDF5_generic( string filename, arraydata<double> *src, int img_rank, hsize_t *dims, int internalType, int verbose ) {
		if ( !src ){
			cerr << "ERROR. In arraydataIO::writeToHDF5_generic. No source data. Could not write to file " << filename << endl;
			return 1;
		}
		
		//output size
		int n = src->size();
		
		hid_t dataspaceID = H5Screate_simple(img_rank, dims, NULL);
		
		std::ostringstream info;
		info << n << " entries written to " << filename << " ";
				
		hid_t memTypeID;				//datatype in memory (not necessarily in the file)
		hid_t memSpaceID = H5S_ALL;		//dataspace in memory, H5S_ALL is default --> the whole dataspace in memory will be used for I/O
		hid_t fileSpaceID = H5S_ALL;	//dataspace in file, H5S_ALL is default --> the whole dataspace in the fill is selected for I/O
		hid_t xfer_prp = H5P_DEFAULT;	//transfer property list, H5P_DEFAULT is default
		void *data = NULL;				//data buffer
		
		// translate the internalType variable to HDF5 types
		info << "(data type ";
		switch(internalType){
			case 0:
				memTypeID = H5T_NATIVE_DOUBLE;
				info << "H5T_NATIVE_DOUBLE";
				break;
			case 1:
				memTypeID = H5T_NATIVE_FLOAT;
				info << "H5T_NATIVE_FLOAT";
				break;
			case 2:
				memTypeID = H5T_NATIVE_INT;
				info << "H5T_NATIVE_INT";
				break;
			case 3:
				memTypeID = H5T_STD_I16LE;
				info << "H5T_STD_I16LE";
				break;
			case 4:
				memTypeID = H5T_NATIVE_LONG;
				info << "H5T_NATIVE_LONG";
				break;
			default:
				memTypeID = H5T_NATIVE_DOUBLE;
				info << "H5T_NATIVE_DOUBLE";
		}
		info << ").";
		
	
		// then allocate data buffer to write	
		if (internalType == 0){
			data = new double[n];
		}else if (internalType == 1){
			data = new float[n];
		}else if (internalType == 2){
			data = new int[n];
		}else if (internalType == 3){
			data = new int16_t[n];
		}else if (internalType == 4){
			data = new long[n];		
		}else{ 								//default = doubles
			cerr << "Warning in  arraydataIO::writeToHDF5. Type not known. Continuing with H5T_NATIVE_DOUBLE." << endl;
			memTypeID = H5T_NATIVE_DOUBLE;
			data = new double[n];
		}
		
		// data assignment and output buffer initialization.
		for( int i = 0; i < n; i++){
			if (memTypeID == H5T_NATIVE_DOUBLE) { 
				((double*)data)[i] = double( src->get_atIndex(i) );
			} else if (memTypeID == H5T_NATIVE_FLOAT) { 
				((float*)data)[i] = float( src->get_atIndex(i) );
			} else if (memTypeID == H5T_NATIVE_INT) { 
				((int*)data)[i] = int( src->get_atIndex(i) );
			} else if (memTypeID == H5T_STD_I16LE) { 
				((int16_t*)data)[i] = int16_t( src->get_atIndex(i) );
			} else if (memTypeID == H5T_NATIVE_LONG) { 
				((long*)data)[i] = long( src->get_atIndex(i) );
			} else {
				((double*)data)[i] = src->get_atIndex(i); 
			}
		}

//		for( int j = 0; j < ny; j++){
//			for( int i = 0; i < nx; i++){
//				if (memTypeID == H5T_NATIVE_DOUBLE) { 
//					((double*)data)[j*nx+i] = double( src->get(i, j) );
//				} else if (memTypeID == H5T_NATIVE_FLOAT) { 
//					((float*)data)[j*nx+i] = float( src->get(i, j) );
//				} else if (memTypeID == H5T_NATIVE_INT) { 
//					((int*)data)[j*nx+i] = int( src->get(i, j) );
//				} else if (memTypeID == H5T_STD_I16LE) { 
//					((int16_t*)data)[j*nx+i] = int16_t( src->get(i, j) );
//				} else if (memTypeID == H5T_NATIVE_LONG) { 
//					((long*)data)[j*nx+i] = long( src->get(i, j) );
//				} else {
//					((double*)data)[j*nx+i] = src->get(i, j); 
//				}
//			}
//		}
		
		// returns 'file' handler
		hid_t fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, xfer_prp, xfer_prp);
		
		// create group "data"
		string groupname = "data";
		hid_t groupID = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, xfer_prp, xfer_prp);
		if (groupID < 0){
			std::cerr << "Error in writeToHDF5. Couldn't create group " << groupname << endl;
			H5Fclose(fileID);
			return 1;
		}
		
		//create dataset in dataspace
		std::string dataset_name = groupname+"/data";
		hid_t datasetID = H5Dcreate1(fileID, dataset_name.c_str(), memTypeID, dataspaceID, xfer_prp);
		if (datasetID < 0) {
			cerr << "Error in arraydataIO::writeToHDF5. Couldn't create dataset." << endl;
			H5Fclose(fileID);
			return 1;
		}
		
		// write data to the dataset
		herr_t status = H5Dwrite(datasetID, memTypeID, memSpaceID, fileSpaceID, xfer_prp, data);
		if (status < 0){ 
			cerr << "Error in arraydataIO::writeToHDF5. Could not write HDF5. " << endl;
			return 2;
		}		
		
		//close all open structures 
		H5Dclose(datasetID);
		H5Sclose(dataspaceID);
		H5Gclose(groupID);
		H5Fclose(fileID);

		if (memTypeID == H5T_NATIVE_DOUBLE) { 
			delete[] ((double*)data);
		} else if (memTypeID == H5T_NATIVE_FLOAT) { 
			delete[] ((float*)data);
		} else if (memTypeID == H5T_NATIVE_INT) { 
			delete[] ((int*)data);
		} else if (memTypeID == H5T_STD_I16LE) { 
			delete[] ((int16_t*)data);
		} else if (memTypeID == H5T_NATIVE_LONG) { 
			delete[] ((long*)data);
		} else { 
			cerr << "Error in arraydataIO::writeToHDF5. Unsupported data format!" << endl;
			return 1;
		}

		if (verbose){
			cout << info.str() << endl;
		}
		return 0;
	}
			
	//-------------------------------------------------------------- writeToHDF5 (1D case)
	int arraydataIO::writeToHDF5( string filename, array1D<double> *src, int internalType ) const {
		if ( !src ){
			cerr << "ERROR. In arraydataIO::writeToHDF5 (1D). No source data. Could not write to file " << filename << endl;
			return 1;
		}
		
		//create dataspace of the right dimensions, pass that down to the more generic function
		int n = src->dim1();
		int img_rank = 1;
		hsize_t	dims[img_rank];
		dims[0] = n;
		
		int fail = writeToHDF5_generic( filename, src, img_rank, dims, internalType, verbose() );

		return fail;
	}

	//-------------------------------------------------------------- writeToHDF5 (2D case)
	int arraydataIO::writeToHDF5( string filename, array2D<double> *src, int internalType ) const {
		if ( !src ){
			cerr << "ERROR. In arraydataIO::writeToHDF5 (2D). No source data. Could not write to file " << filename << endl;
			return 1;
		}
		
		if (transposeForIO()) { src->transpose(); }

		//create dataspace of the right dimensions, pass that down to the more generic function
		int nx = src->dim1();
		int ny = src->dim2();
		int img_rank = 2;
		hsize_t	dims[img_rank];
		dims[0] = ny;
		dims[1] = nx;		

		int fail = writeToHDF5_generic( filename, src, img_rank, dims, internalType, verbose() );
		
		if (transposeForIO()) { src->transpose(); }		//transpose back before returning (preserve original data)
		
//		if (debug) {
//			cout << "nx: " << nx << ", ny: " << ny << endl; 
//			
//			//observe original data
//			cout << "original data: " << endl;
//			for( int j = 0; j < ny; j++){
//				cout << "j" << j << ": ";
//				for( int i = 0; i < nx; i++){
//					cout << src->get(i, j) << " ";
//				}
//				cout << endl;
//			}
//
//			//observe data stream to be written
//			cout << "data buffer:" << endl;
//			for( int j = 0; j < ny; j++){
//				cout << "j" << j << ": ";
//				for( int i = 0; i < nx; i++){
//					if (memTypeID == H5T_NATIVE_DOUBLE) { 
//						cout << ((double*)data)[j*nx+i] << " ";
//					} else if (memTypeID == H5T_NATIVE_FLOAT) { 
//						cout << ((float*)data)[j*nx+i] << " ";
//					} else if (memTypeID == H5T_NATIVE_INT) { 
//						cout << ((int*)data)[j*nx+i] << " ";
//					} else if (memTypeID == H5T_STD_I16LE) { 
//						cout << ((int16_t*)data)[j*nx+i] << " ";
//					} else if (memTypeID == H5T_NATIVE_LONG) { 
//						cout << ((long*)data)[j*nx+i] << " ";
//					} else { 
//						cout << "Error in arraydataIO::writeToHDF5. Unsupported data format!" << endl;
//						return 1;
//					}
//				}//for i
//				cout << endl;
//			}//for j
//		}//debug
		return fail;
	}
	
#else //define empty dummy functions
	int arraydataIO::readFromHDF5( string filename, array2D<double> *&dest ) const { 
		cerr << "======== WARNING! Dummy function. Nothing read from HDF5. ======== " << endl;
		return 1;
	}
	int arraydataIO::writeToHDF5( string filename, array2D<double> *src, int type, int debug ) const { 
		cerr << "======== WARNING! Dummy function. Nothing written to HDF5. ======== " << endl; 
		return 1;
	}
#endif





//****************************************************************************************************************
//****************************************************************************************************************
//*** ASCII                                                                                                    ***
//****************************************************************************************************************
//****************************************************************************************************************


//-------------------------------------------------------------- readFromASCII (1D)
//
//-------------------------------------------------------------- 
int arraydataIO::readFromASCII( std::string filename, array1D<double> *&dest ) const {
	array2D<double> *twoD = 0;
	int fail = readFromASCII( filename, twoD );
	if (!fail){
		delete dest;
		dest = new array1D<double>( twoD );
	}else{
		return fail;
	}
	delete twoD;
	return 0;
}

//-------------------------------------------------------------- readFromASCII (12D)
//
//-------------------------------------------------------------- 
int arraydataIO::readFromASCII( std::string filename, array2D<double> *&dest ) const {
	if (verbose()) {
		cout << "Reading from ASCII file " << filename << endl;
	}
	std::ifstream fin( filename.c_str() );
	if (!fin.fail()) {
		unsigned int linecount = 0; 
		unsigned int colcount = 0;

		std::vector< std::vector<string> > matrix_str;
		while (!fin.eof()) {
			std::vector<string> numbers_str;
			char cbuffer[4096];
			fin.getline(cbuffer, 4096, '\n');
			std::istringstream iss(cbuffer);
			copy(std::istream_iterator<string>(iss), 
					std::istream_iterator<string>(), 
					std::back_inserter<std::vector<string> >(numbers_str) );
			matrix_str.push_back( numbers_str );
			if ( numbers_str.size() > colcount ){
				colcount = (unsigned int) numbers_str.size(); 
			}
			linecount++;
		}//while
		linecount--;
		
		if (verbose()>1) {
			cout << "read lines: " << linecount << ", columns: " << colcount << endl;
		}
		
		// transfer to 'dest'
		delete dest;
		dest = new array2D<double>( linecount, colcount );
		for( unsigned int row = 0; row < linecount; row++ ){
			for( unsigned col = 0; col < colcount; col++ ){
				std::istringstream isst( matrix_str.at(row).at(col) );
				double val = 0.;
				isst >> val;
				dest->set( row, col, val );
			}	
		}
		
		fin.close();
	}else{
		cerr << "ERROR in arraydataIO::readFromASCII. Couldn't read from file " << filename << ". " << endl;
		return 1;
	}
	return 0;
}

//-------------------------------------------------------------- writeToASCII (1D)
//
//-------------------------------------------------------------- 
int arraydataIO::writeToASCII( std::string filename, array1D<double> *src, int format ) const {
	if ( !src ){
		cerr << "ERROR in arraydataIO::writeToASCII (1D). No source data. Could not write to file " << filename << endl;
		return 1;
	}
	
	if (src->size() == 0) {
		cerr << "ERROR in arraydataIO::writeToASCII (1D). Array size is zero." << endl;
		return 2;
	}

	std::ofstream fout( filename.c_str() );
	switch (format){
		case 1:
			fout << src->getASCIIdataAsColumn();
			break;
		case 2:
			fout << src->getASCIIdataAsRow();
			break;
		default:
			fout << src->getASCIIdata(0);
	}
	fout.close();
	return 0;
}

//-------------------------------------------------------------- writeToASCII (2D)
//
//-------------------------------------------------------------- 
int arraydataIO::writeToASCII( std::string filename, array2D<double> *src, int format ) const {
	if ( !src ){
		cerr << "ERROR in arraydataIO::writeToASCII (2D). No source data. Could not write to file " << filename << endl;
		return 1;
	}
	
	if (src->size() == 0) {
		cerr << "ERROR in arraydataIO::writeToASCII (2D). Array size is zero." << endl;
		return 2;
	}
	

	std::ofstream fout( filename.c_str() );
	switch (format){
		case 1:
			fout << src->getASCIIdataAsColumn();
			break;
		case 2:
			fout << src->getASCIIdataAsRow();
			break;
		default:
			fout << src->getASCIIdata(0);
	}
	fout.close();
	return 0;
}

//-------------------------------------------------------------- writeToASCII (3D)
//
//-------------------------------------------------------------- 
int arraydataIO::writeToASCII( std::string filename, array3D<double> *src, int format ) const {
	if ( !src ){
		cerr << "ERROR in arraydataIO::writeToASCII (3D). No source data. Could not write to file " << filename << endl;
		return 1;
	}
	
	if (src->size() == 0) {
		cerr << "ERROR in arraydataIO::writeToASCII (3D). Array size is zero." << endl;
		return 2;
	}	


	std::ofstream fout( filename.c_str() );
	switch (format){
		case 1:
			fout << src->getASCIIdataAsColumn();
			break;
		case 2:
			fout << src->getASCIIdataAsRow();
			break;
		default:
			fout << src->getASCIIdata(0);
	}
	fout.close();
	return 0;
}






//****************************************************************************************************************
//****************************************************************************************************************
//*** general functionality                                                                                    ***
//****************************************************************************************************************
//****************************************************************************************************************

bool arraydataIO::transposeForIO() const{
	return p_transposeForIO;
}

void arraydataIO::setTransposeForIO( bool t ){
	p_transposeForIO = t;
}

void arraydataIO::setVerbose( int level ){
	p_verbose = level;
}

int arraydataIO::verbose() const{
	return p_verbose;
}



//-----------------------------------------------------test pattern
void arraydataIO::generateTestPattern( array2D<double> *img, int type ){
	if (img == 0){
		cerr << "Error in arraydataIO::generateTestPattern. Image not allocated." << endl;
	}
	
	if (img->dim1()==0 || img->dim2()==0)
		cout << "WARNING in generateTestPattern. One or more dimensions are zero." << endl;


	cout << "Writing test pattern type " << type << ": ";
	switch (type) {
		case 0:											
			{   
				cout << "2D sinusoidal ";
				double amplitude = 20000;
				double periodX = img->dim1()/1;
				double periodY = img->dim2()/2;
				for (int i = 0; i < img->dim1(); i++) {
					for (int j = 0; j < img->dim2(); j++) {
						img->set(i, j, amplitude/2*(sin(2*M_PI*i/periodX)*cos(2*M_PI*j/periodY)+1) );
					}
				}
			}
			break;
		case 1:											
			{
				cout << "increment by absolute array index ";
				double val = 0;
				for (int i = 0; i < img->size(); i++) {
					if (val >= 65535) {
						val = 0; 
					}
					img->set_atIndex(i, val);
					val += 1;
				}
			}
			break;
		case 2:											
			{
				cout << "2D centro-symmetric sine ";
				double amplitude = 20000;
				double period = img->dim1()/5;
				double centerX = img->dim1()/2;
				double centerY = img->dim2()/2;
				cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < img->dim1(); i++) {
					for (int j = 0; j < img->dim2(); j++) {
						double x = i - centerX;
						double y = j - centerY;
						double r = sqrt( x*x + y*y );
						img->set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 ) );
					}
				}
			}
			break;
		case 3:											
			{
				cout << "2D centro-symmetric sine with circular modulation ";
				double amplitude = 20000;
				double period = img->dim1()/5;
				double periodMod = img->dim1()/10;
				double centerX = img->dim1()/2;
				double centerY = img->dim2()/2;
				cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < img->dim1(); i++) {
					for (int j = 0; j < img->dim2(); j++) {
						double x = i - centerX;
						double y = j - centerY;
						double r = sqrt( x*x + y*y );
						double phi = atan(y/x);
						img->set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 + 1/10*(tan(2*M_PI*phi/periodMod)+1) ) );
					}
				}
			}
			break;
		case 4:											
			{
				cout << "2D centro-symmetric sine with straight modulation ";
				double amplitude = 20000;
				double period = img->dim1()/5;
				double periodMod = img->dim1()/10;
				double centerX = img->dim1()/2;
				double centerY = img->dim2()/2;
				cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < img->dim1(); i++) {
					for (int j = 0; j < img->dim2(); j++) {
						double x = i - centerX;
						double y = j - centerY;
						double r = sqrt( x*x + y*y );
						img->set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 + cos(2*M_PI*i/periodMod)+1 ) );
					}
				}
			}
			break;
		default:                                        
			cout << "zeros";
			img->zeros();
			break;
	}
	cout << endl;
}





