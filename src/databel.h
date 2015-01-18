/*
 * Copyright (c) 2010-2015, Diego Fabregat-Traver and Paolo Bientinesi.
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

#ifndef DATABEL_H
#define DATABEL_H

enum datatype { UNSIGNED_SHORT_INT_TYPE = 1,
                SHORT_INT_TYPE,
                UNSIGNED_INT_TYPE,
                INT_TYPE,
                FLOAT_TYPE,
                DOUBLE_TYPE,
                SIGNED_CHAR_TYPE,
                UNSIGNED_CHAR_TYPE };

static char *type_string[8] = {
	"UNSIGNED SHORT INT",
	"SHORT_INT",
	"UNSIGNED_INT",
	"INT",
	"FLOAT",
	"DOUBLE",
	"SIGNED CHAR",
	"UNSIGNED CHAR"
};

#define TYPE2STR(TYPE) (type_string[TYPE-1])

#define NAMELENGTH 32
#define RESERVEDSPACE 5

typedef struct databel_fvi_header
{
	unsigned short int type;
	unsigned int nelements;
	unsigned int numObservations;
	unsigned int numVariables;
	unsigned int bytesPerRecord;
	unsigned int bitsPerRecord;
	unsigned int namelength;
	unsigned int reserved[RESERVEDSPACE];
} databel_fvi_header;

typedef struct databel_fvi
{
	databel_fvi_header  fvi_header;
	char               *fvi_data;
} databel_fvi;

struct databel_fvi * load_databel_fvi( char *path );
void free_databel_fvi( struct databel_fvi **fvi );

#endif // DATABEL_H
