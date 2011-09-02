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

#include <cmath>


arraydataIO::arraydataIO(){
	p_transpose = true;			//should be 'true', generally. see note in header
}

arraydataIO::~arraydataIO(){

}

#ifdef ARRAYDATAIO_EDF
	#include "edf.h"				// needed to read/write EDF files (ESRF Data Format)
	
	//-------------------------------------------------------------- readFromEDF
	// implementation borrowed following matrixdata::in(string fn, datatype_t &datatype) 
	// of the tomo/matlib package (see http://xray-lens.de )
	//--------------------------------------------------------------------------
	int arraydataIO::readFromEDF( string filename, array2D *&dest ) const{
		int retval = 0;
		ns_edf::edf *file = new ns_edf::edf;
		if (dest->data())
		{
			file->read_header(filename);
			
			int dim1 = (int)file->get_Dim_x();
			int dim2 = (int)file->get_Dim_y();
			
			cout << "Reading " << filename << " (EDF file type " << file->get_FileType_str() << ")"
				<< ", dimensions (" << dim1 << ", " << dim2 << ")" << endl;
			
			double *temp = new double[dim1*dim2];
			
			int fail = file->read_data( temp, filename );
			if ( fail ){
				cout << "ERROR. In arraydataIO::readFromEDF - Could not load EDF data file '" << filename << "'." << endl;
				return 1; // Could not open file fn   (by JP)
			}
			
			//feed back into dest
			delete dest;
			dest = new array2D( dim1, dim2 );
			dest->arraydata::copy( temp, dim1*dim2);
			
			delete []temp;
			
			if (transpose()){
				dest->transpose();
			}
		}
		else
		{
			cerr << "ERROR. arraydataIO::readFromEDF - dest data not allocated." << endl;
			retval = 1;
		}
		delete file;
		return retval; 
	}

	//-------------------------------------------------------------- writeToEDF
	//
	//--------------------------------------------------------------
	int arraydataIO::writeToEDF( string filename, array2D *src ) const {
		int retval = 0;
		bool flipByteOrder = false;
		
		// write scaled data?
		// SF_UNSCALED = 0:  binary values written represent actual data values (a range from 0 to 65535 is possible)
		// SF_SCALED	= 1: scaling max and min are written to header, data values always span full range from 0 to 65535 (2^16-1)
		// SF_RETAIN	= 2: if the scaling factor was read from an input file, keep it what it was
		ns_edf::scaled_t scaleOut = ns_edf::SF_SCALED;				
		
		if ( src->data() ){
			if (transpose()){//transpose before writing
				src->transpose();
			}
			ns_edf::edf *file = new ns_edf::edf;
			
			file->set_Dim_x( src->dim1() );
			file->set_Dim_y( src->dim2() );
			file->set_ScalingMin( src->calcMin() );
			file->set_ScalingMax( src->calcMax() );
			file->set_FilenameOut( filename );
			
	//		file->set_FileType( FTS_MATRIX_DATA );
	//		file->set_DataType( datatype );
			file->set_ScaledFlag( scaleOut );
			retval = file->write( src->data(), flipByteOrder );
			delete file;
			
			if (transpose()){//transpose back before returning
				src->transpose();
			}
			return 0;
		}else{
			cerr << "ERROR. In matrixdata::write() - Matrix not allocated." << endl;
			return 1;
		}
		return retval;
	}

#else
	int arraydataIO::readFromEDF( string filename, array2D *&dest ) const{
		cout << "======== WARNING! Dummy function. Nothing read from EDF. ======== " << endl; 
		return 1;
	}
	int arraydataIO::writeToEDF( string filename, array2D *src ) const {
		cout << "======== WARNING! Dummy function. Nothing written to EDF.======== " << endl; 
		return 1;
	}
#endif


#ifdef ARRAYDATAIO_TIFF
	#include <tiffio.h>				// needed to read/write tiff files (Tagged Image File Format)
	
	//-------------------------------------------------------------- readFromTiff
	// (needs TIFFLIB installation)
	//
	// implementation borrowed following matrixdata::tiffin 
	// of the tomo/matlib package (see http://xray-lens.de )
	//--------------------------------------------------------------------------
	int arraydataIO::readFromTiff( string filename, array2D *&dest ) const {
		int retval = 0;
		TIFF *tiff= TIFFOpen( filename.c_str(), "r" );

		if(tiff){
			uint32 width, height;
			int npixels;
			uint32* raster;

			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);

			//allocate new data to write image to		
			delete dest;
			dest = new array2D(width, height);

			cout << "Reading " << filename << " (TIFF file)"
				<< ", dimensions (" << dest->dim1() << ", " << dest->dim2() << ")" << endl;

			npixels = width * height;
			raster = (uint32*) _TIFFmalloc( npixels * sizeof(uint32) );
			if (raster != NULL){
				if ( TIFFReadRGBAImage(tiff, width, height, raster, 0) ){
					dest->setDim1(width);
					dest->setDim2(height);
					for(unsigned int j= 0; j < height; j++){		//y
						for(unsigned int i= 0; i < width; i++){		//x	
//							dest->set( i, j, 1.* (raster[j*width+i] & 0x000000ff));			// max value 255
							dest->set( i, j, 1.* (raster[j*width+i] & 0x0000ffff));			// max value 65535
						}
					}
				}
				_TIFFfree(raster);
			}
			TIFFClose(tiff);
			
			if (transpose()){
				dest->transpose();
			}
			
			dest->flipud();
		} else {
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
	int arraydataIO::writeToTiff( string filename, array2D *src, int scaleFlag, int verbose ) const {

		if (src->size() == 0) {
			cerr << "Error in writeToTiff! Array size is zero." << endl;
		}
		
		if (!src->data()) {
			cerr << "Error in writeToTiff! No data in array." << endl;
		}
		
		TIFF *out = TIFFOpen(filename.c_str() ,"w");
		if(out)
		{
			if (transpose()){
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
			for (uint32 i = 0; i < height; i++)
			{
				for (unsigned int j = 0; j < width; j++)
				{
					if (scaleFlag) {										//scale image to full range
						tifdata[j] = (uint16) floor( 65535.* (src->get(j,i)-MinValue) / MaxMinRange );
					} 
					else {										
						if( src->get(j,i) < 0 ) {
							tifdata[j] = 0;									// cut off values smaller than 0
						}
						else if ( src->get(j,i) > 65535 ) {
							tifdata[j] = 65535;								// ...and values larger than 2^16-1
						}
						else{ 
							tifdata[j] = (uint16) floor(src->get(j,i));		// ... force everything else to be integer
						}
					}
				}
				
				TIFFWriteScanline(out, tifdata, i, 0);
			}
			
			
			if (verbose){
				cout << "Tiff image '" << filename << "' written to disc." << endl;
				if (scaleFlag) cout << "Scaled output. "; else cout << "Unscaled output. ";
				cout << "Data min: " << MinValue << ", max: " << MaxValue << ", range: " << MaxMinRange << endl;
			}
			
			delete[] tifdata;
			TIFFClose(out);
			
			if (transpose()){//transpose back before returning
				src->transpose();
			}
			return 0;
		}else{
			cerr << "Error in array2D::writeToTiff! Could not open '" << filename << "' for writing!" << endl;
			return 1;	
		}
	}

#else
	int arraydataIO::readFromTiff( string filename, array2D *&dest ) const {
		cout << "======== WARNING! Dummy function. Nothing read from Tiff. ======== " << endl; 
		return 1;
	}
	int arraydataIO::writeToTiff( string filename, array2D *src, int scaleFlag, int verbose ) const {
		cout << "======== WARNING! Dummy function. Nothing written to Tiff. ======== " << endl; 
		return 1;
	}
#endif



#ifdef ARRAYDATAIO_HDF5
	#include <hdf5.h>				// needed to read/write HDF5 files (Tagged Image File Format)
	
	//-------------------------------------------------------------- readFromHDF5
	//
	//-------------------------------------------------------------- 
	int arraydataIO::readFromHDF5( string filename, array2D *&dest ) const {
		cout << "Reading " << filename << " (HDF5 file)" << endl;
		
		string dataset_name = "data/data";

		
		// open file and dataset using the default properties
		hid_t fileID = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		hid_t datasetID = H5Dopen (fileID, dataset_name.c_str(), H5P_DEFAULT);
		
		hid_t typeID = H5Dget_type(datasetID);
		H5T_class_t t_class = H5Tget_class(typeID);

		H5T_order_t order = H5Tget_order(typeID);
		if (order == H5T_ORDER_LE){
			cout << "Little endian order" << endl;
		}else if (order == H5T_ORDER_BE){
			cout << "Big endian order" << endl;
		}

		size_t size = H5Tget_size(typeID);

		hid_t dataspace = H5Dget_space(datasetID);    /* dataspace handle */
		int rank = H5Sget_simple_extent_ndims(dataspace);
		hsize_t dims_out[2];
		int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
		int ny = (int)(dims_out[0]);
		int nx = (int)(dims_out[1]);
		cout << "rank " << rank << ", dimensions nx:" << nx << " x ny:" << ny << ", status " << status_n << endl;

		//allocate data to hold what comes out of the file
		void *buffer = NULL;
		int type = 0;
		if (t_class == H5T_FLOAT){
			cout << "Data set has FLOAT type, data size is " << (int)size;
			if (size == sizeof(float)){
				cout << " --> float" << endl;
				buffer = new float[nx*ny];
				type = 1;
			}else if (size == sizeof(double)){
				cout << " --> double" << endl;
				buffer = new double[nx*ny];
				type = 0;
			}
		}else if (t_class == H5T_INTEGER){
			cout << "Data set has INTEGER type, data size is " << (int)size;
			type = 2;
			if (size == sizeof(int)){
				cout << " --> int" << endl;
				buffer = new int[nx*ny];
			}else if (size == sizeof(short int)){
				cout << " --> short" << endl;
				buffer = new short int[nx*ny];
			}else if (size == sizeof(long int)){
				cout << " --> long" << endl;
				buffer = new long int[nx*ny];
			}
		}else{
			cerr << "Error in arraydataIO::readFromHDF5. Unknown HDF5 class type '" << t_class << "'" << endl;
		}

		// read the data using the default properties.
	//    hid_t status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);
		hid_t status = H5Dread (datasetID, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
		
		// close and release resources.
		status = H5Dclose(datasetID);
		status = H5Fclose(fileID);
		
		//copy data back to dest
		delete dest;
		dest = new array2D( nx, ny );


		//for some reason, the dereferencing of anything else than 'int' seems to fail
		//need to fix this as this throws away all floating point accurracy
		dest->arraydata::copy( (int*)buffer, nx*ny);

		if (type == 0) { 
			//dest->arraydata::copy( (double*)buffer, nx*ny);
			delete[] ((double*)buffer);
		} else if (type == 1) { 
			//dest->arraydata::copy( (float*)buffer, nx*ny);
			delete[] ((float*)buffer);
		} else if (type == 2) { 
			//dest->arraydata::copy( (int*)buffer, nx*ny);
			delete[] ((int*)buffer);
		} else { 
			cerr << "Error in arraydataIO::readFromHDF5. Type not found. " << endl;
		}
		
		if (transpose()){
			dest->transpose();
		}
		
		return 0;
	}

	//-------------------------------------------------------------- writeToHDF5
	// dataType = 0 --> write as doubles (default) 
	// dataType = 1 --> write as float
	// dataType = 2 --> write as int
	//-------------------------------------------------------------- 
	int arraydataIO::writeToHDF5( string filename, array2D *src, int type, int debug ) const {

		if (transpose()){
			src->transpose();
		}
		
		int nx = src->dim1();
		int ny = src->dim2();

		hid_t typeID;
		void *data;
		
		// allocate data buffer to write	
		if (type == 0){
			typeID = H5T_NATIVE_DOUBLE;
			data = new double[nx*ny];
		}else if (type == 1){
			typeID = H5T_NATIVE_FLOAT;
			data = new float[nx*ny];
		}else if (type == 2){
			typeID = H5T_NATIVE_INT;	//H5T_STD_I16LE creates weird results --> look into this later (maybe)....
			data = new int[nx*ny];
		}else{ 								//default = doubles
			cerr << "Warning in  arraydataIO::writeToHDF5. Type not known. Continuing with H5T_NATIVE_DOUBLE." << endl;
			typeID = H5T_NATIVE_DOUBLE;
			data = new double[nx*ny];
		}
		
		// data assignment and output buffer initialization.
		for( int j = 0; j < ny; j++){
			for( int i = 0; i < nx; i++){
				if (typeID == H5T_NATIVE_DOUBLE) { 
					((double*)data)[j*nx+i] = double( src->get(i, j) );
				} else if (typeID == H5T_NATIVE_FLOAT) { 
					((float*)data)[j*nx+i] = float( src->get(i, j) );
				} else if (typeID == H5T_NATIVE_INT) { 
					((int*)data)[j*nx+i] = int( src->get(i, j) );
				} else { 
					((double*)data)[j*nx+i] = src->get(i, j); 
				}
			}
		}
		
		// create file, ANTON_CHECK
		// returns 'file' handler
		hid_t fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		
		
		// create group "data"
		string groupname = "data";
		hid_t groupID = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (groupID < 0){
			std::cerr << "Error in writeToHDF5. Couldn't create group " << groupname << endl;
			H5Fclose(fileID);
			return 1;
		}

		
		//create dataspace
		int img_rank = 2;
		hsize_t	dims[img_rank];
		dims[0] = ny;
		dims[1] = nx;
		hid_t dataspaceID = H5Screate_simple(img_rank, dims, NULL);
		
		
		//create dataset in dataspace
		std::string dataset_name = groupname+"/data";
		hid_t datasetID = H5Dcreate1(fileID, dataset_name.c_str(), typeID, dataspaceID, H5P_DEFAULT);
		if (datasetID < 0) {
			cerr << "Error in arraydataIO::writeToHDF5. Couldn't create dataset." << endl;
			H5Fclose(fileID);
		}
		
		
		// write data to the dataset
		herr_t status = H5Dwrite(datasetID, typeID, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		if (status < 0){ 
			cerr << "Error in arraydataIO::writeToHDF5. Could not write HDF5. " << endl;
		}
		
		
		//close all open structures 
		H5Dclose(datasetID);
		H5Sclose(dataspaceID);
		
		H5Gclose(groupID);
		H5Fclose(fileID);
		
		if (debug) {
			cout << "nx: " << nx << ", ny: " << ny << endl; 
			
			//observe original data
			cout << "original data: " << endl;
			for( int j = 0; j < ny; j++){
				cout << "j" << j << ": ";
				for( int i = 0; i < nx; i++){
					cout << src->get(i, j) << " ";
				}
				cout << endl;
			}

			//observe new data stream to be written
			cout << "data buffer:" << endl;
			for( int j = 0; j < ny; j++){
				cout << "j" << j << ": ";
				for( int i = 0; i < nx; i++){
					if (typeID == H5T_NATIVE_DOUBLE) { 
						cout << ((double*)data)[j*nx+i] << " ";
					} else if (typeID == H5T_NATIVE_FLOAT) { 
						cout << ((float*)data)[j*nx+i] << " ";
					} else if (typeID == H5T_NATIVE_INT) { 
						cout << ((int*)data)[j*nx+i] << " ";
					} else { 
						cout << ((double*)data)[j*nx+i] << " ";
					}
				}//for i
				cout << endl;
			}//for j
		}//debug

		if (typeID == H5T_NATIVE_DOUBLE) { 
			delete[] ((double*)data);
		} else if (typeID == H5T_NATIVE_FLOAT) { 
			delete[] ((float*)data);
		} else if (typeID == H5T_NATIVE_INT) { 
			delete[] ((int*)data);
		} else { 
			delete[] ((double*)data);
		}
		
		if (transpose()){//transpose back before returning
			src->transpose();
		}
		return 0;
	}
#else //define empty dummy functions
	int arraydataIO::readFromHDF5( string filename, array2D *&dest ) const { 
		cout << "======== WARNING! Dummy function. Nothing read from HDF5. ======== " << endl;
		return 1;
	}
	int arraydataIO::writeToHDF5( string filename, array2D *src, int type, int debug ) const { 
		cout << "======== WARNING! Dummy function. Nothing written to HDF5. ======== " << endl; 
		return 1;
	}
#endif



bool arraydataIO::transpose() const{
	return p_transpose;
}

void arraydataIO::setTranspose( bool t ){
	p_transpose = t;
}




