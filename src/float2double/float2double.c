/*
 * Copyright (c) 2010-2013, Diego Fabregat-Traver and Paolo Bientinesi.
 * All rights reserved.
 *
 * This file is part of OmicABEL.
 * 
 * OmicABEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * OmicABEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with OmicABEL. If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * Coded by:
 *   Diego Fabregat-Traver (fabregat@aices.rwth-aachen.de)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../wrappers.h"
#include "../databel.h"

#define MB (1L<<20)
#define STR_BUFFER_SIZE 256

int main( int argc, char *argv[] )
{
	char  fin_path_fvi[STR_BUFFER_SIZE], 
	      fin_path_fvd[STR_BUFFER_SIZE], 
		 fout_path_fvi[STR_BUFFER_SIZE],
		 fout_path_fvd[STR_BUFFER_SIZE];
	FILE *fin, *fout;
	struct databel_fvi *databel_in, *databel_out;

	float *datain;
	double *dataout;
	size_t buff_size = 256*MB;

	long long int nelems;
	int nelems_in_buff, nelems_to_write;
	int header_data_size;

	int i, j, out;

	if ( argc != 3 )
	{
		fprintf( stderr, "Usage: %s floatFileIn doubleFileOut\n", argv[0] );
		exit( EXIT_FAILURE );
	}

	snprintf(  fin_path_fvi, STR_BUFFER_SIZE, "%s.fvi", argv[1] ); 
	snprintf(  fin_path_fvd, STR_BUFFER_SIZE, "%s.fvd", argv[1] ); 
	snprintf( fout_path_fvi, STR_BUFFER_SIZE, "%s.fvi", argv[2] ); 
	snprintf( fout_path_fvd, STR_BUFFER_SIZE, "%s.fvd", argv[2] ); 

	// FVI files
	databel_in = load_databel_fvi( fin_path_fvi );
	if ( databel_in->fvi_header.type != FLOAT_TYPE )
	{
		fprintf( stderr, "Input databel file(s) %s should include \"float\" data\n", argv[1]);
		exit( EXIT_FAILURE );
	}
	databel_out = (databel_fvi *) fgls_malloc( sizeof(databel_fvi) );
	// Header
	databel_out->fvi_header.type = DOUBLE_TYPE;
	databel_out->fvi_header.nelements       = databel_in->fvi_header.nelements;
	databel_out->fvi_header.numObservations = databel_in->fvi_header.numObservations;
	databel_out->fvi_header.numVariables    = databel_in->fvi_header.numVariables;
	databel_out->fvi_header.bytesPerRecord  = sizeof( double );
	databel_out->fvi_header.bitsPerRecord   = databel_out->fvi_header.bytesPerRecord * 8;
	databel_out->fvi_header.namelength      = databel_in->fvi_header.namelength;
	for ( i = 0; i < RESERVEDSPACE; i++ )
		databel_out->fvi_header.reserved[i] = '\0';
	// Labels
	header_data_size = (databel_out->fvi_header.numVariables + databel_out->fvi_header.numObservations ) * 
	                    databel_out->fvi_header.namelength * sizeof(char);
	databel_out->fvi_data = (char *) fgls_malloc ( header_data_size );
	memcpy( databel_out->fvi_data, databel_in->fvi_data, header_data_size );

	// Write
	fout = fgls_fopen( fout_path_fvi, "wb" );
    out = fwrite( &databel_out->fvi_header, sizeof(databel_fvi_header), 1, fout);
    if ( out != 1 ) 
    {
        fprintf(stderr, "Error writing fvi header\n" );
		exit( EXIT_FAILURE );
    }
    out = fwrite( databel_out->fvi_data, 
	              databel_out->fvi_header.namelength * sizeof(char), 
				  databel_out->fvi_header.numVariables + databel_out->fvi_header.numObservations,
				  fout);
    if ( out != (databel_out->fvi_header.numVariables + databel_out->fvi_header.numObservations) ) 
    {
        fprintf(stderr, "Error writing fvi data\n" );
		exit( EXIT_FAILURE );
	}
	fclose( fout );

	// FVD
	fin  = fgls_fopen(  fin_path_fvd, "rb" );
	fout = fgls_fopen( fout_path_fvd, "wb" );
	// buff_size determines the size of the buffer for the "double" array.
	// For the same amount of elements, float needs half the memory space
	datain  = (float *)  fgls_malloc( buff_size / 2 );
	dataout = (double *) fgls_malloc( buff_size );

	nelems = databel_out->fvi_header.numVariables * databel_out->fvi_header.numObservations; // total elems in file
	nelems_in_buff = buff_size / sizeof(double);
	for ( i = 0; i < nelems; i += nelems_in_buff )
	{
		nelems_to_write = ((nelems - i) >= nelems_in_buff) ? nelems_in_buff : nelems - i;
		if ( fread( datain, sizeof(float), nelems_to_write, fin ) != nelems_to_write )
		{
			fprintf( stderr, "Error reading data from %s\n", fin_path_fvd );
			exit( EXIT_FAILURE );
		}
		for ( j = 0; j < nelems_to_write; j++ )
		{
#if 0
			if (isnan( datain[j] ))
				printf("Nan in input:  %lld\n", datain[j]);
#endif
			dataout[j] = (double)datain[j];
#if 0
			if (isnan( dataout[j] ))
				printf("Nan in output: %lld\n", dataout[j]);
#endif
		}
		if ( fwrite( dataout, sizeof(double), nelems_to_write, fout ) != nelems_to_write )
		{
			fprintf( stderr, "Error writing data to %s\n", fout_path_fvd );
			exit( EXIT_FAILURE );
		}
	}
	fclose( fin );
	fclose( fout );
	free( datain );
	free( dataout );
	free_databel_fvi( &databel_in );
	free_databel_fvi( &databel_out );
	
	return 0;
}
