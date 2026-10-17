/*=======================================================================
 * All files in the distribution of the DPS system are Copyright
 * 1996 by the Computational Biology group in the Department of Biological
 * Sciences at Purdue University.  All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that this entire copyright notice is duplicated in all such
 * copies, and that any documentation, announcements, and other materials
 * related to such distribution and use acknowledge that the software was
 * developed by the Computational Biology group in the Department of
 * Biological Sciences at Purdue University, W. Lafayette, IN by Ingo
 * Steller and Michael G. Rossmann. No charge may be made for copies,
 * derivations, or distributions of this material without the express
 * written consent of the copyright holder.  Neither the name of the
 * University nor the names of the authors may be used to endorse or
 * promote products derived from this material without specific prior
 * written permission.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE.
 *======================================================================*/

/*=====================================================================*
 *                                                                     *
 *                         Data Processing Suit                        *
 *                                                                     *
 *                            Utility-Library                          *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                           File: read_frame.c                        *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#if !defined vms && !defined (__APPLE__) && ! ( defined (__FreeBSD__) && defined (__STDC__) )
#include <malloc.h>
#endif
#include "util.h"
#include "util_def.h"
#include "read_frame.h"

/* Definition of the different parameters for the different filetypes */

FrameInfo staticframeinfo[13] = {
	{0, 0, 0, 0, 0, 0, 0},
	{2048, 2560, 0, 9821, FUJI, 0, 9910},
	{0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0},
	{3072, 3072, 0, 65534, APS3000, 512, 65535},
	{2048, 1900, 0, 65534, RAXIS2_100, 4096, 65535},
	{3000, 3000, 0, 65534, RAXIS4_100, 6000, 65535},
	{1152, 1152, 0, 65534, QUANTUM1, 512, 65535},
	{2304, 2304, 0, 65534, QUANTUM4, 512, 65535},
	{1200, 1200, 0, 65535, MAR180, 2400, 65535},
	{2000, 2000, 0, 65535, MAR300, 4000, 65535},
	{3000, 3000, 0, 65535, MAR345, 6000, 65535}
};


/* read_frame(): On a given Filename the routine checks if the file exists 
	and if it is a valid optical density file. It opens the file,  identifies 
	the type by looking on the magic number and the filesize, extracts 
	certain information from the header if possible and fills the 
	frame_info structure. It returns a pointer to the optical density data
	that is normally memory mapped from the file,  and a return value
	for the status. The meaning of the reutrn value is:
	0: ok,  file read.
	1: File could not be opend (does not exists,  no permission)
	2: Format not recognized
	3: memory allocation error	    		*/

int read_frame(char *filename, FrameInfo *frame_info, od_type **od_pointer)

{
	int i;
	int in_file;
	char magic[LENGTH_MAGIC+1];
	struct stat filestat;
	od_type *value;
	int ret_value;
	od_type max,  min;

	char *smv_header;
	char smv_value[80];
	int lhead, naxis, type;
	int axis[10];

#ifdef READ_FRAME_DEBUG
	printf("READ_FRAME entered...\n");
#endif

	/* Open file */

	in_file = open(filename, O_RDONLY);
	if (in_file == -1) {
		printf("Can't open %s!! Exit/n", filename);
		return(1);
	};

	/* Read the first LENGTH_MAGIC bytes and get the filesize */
	read(in_file, magic, LENGTH_MAGIC);
	magic[LENGTH_MAGIC]='\0';
	lseek(in_file, 0, SEEK_SET);
	fstat(in_file, &filestat);
#ifdef READ_FRAME_DEBUG
	printf("MAGIC: [%s] length %d\n", magic, filestat.st_size);
#endif

	/* Try to recognize the filetype by analysing the magic number 
       and the filesize */

	/* FUJI 10bit format */
	if(filestat.st_size == 10485760) {
		printf("Recognized FUJI 10bit...\n");
		/* Copy the frame information */
		*frame_info = staticframeinfo[FUJI];

		ret_value = 0;
	}
	/* APS format (3x3 CCD detector 3072 x 3072 pixels) no support for
       analysis of header until now */
	else if((strcmp("{\nHEADER_B", magic) == 0) && (filestat.st_size == 18874880)) {
		printf("Recognized APS3000...\n");
		/* Copy the frame information */
		*frame_info = staticframeinfo[APS3000];

		ret_value = 0;
	}
	/* RAXIS2 100micron */
	else if((strcmp("R-AXIS2   ", magic) == 0)  && (filestat.st_size == 7786496)) {
		printf("Recognized R-AXIS2 100micron...\n");
		/* Copy the frame information */
		*frame_info = staticframeinfo[RAXIS2_100];

		ret_value = 0;
	}
	/* RAXIS4 100micron */
	else if((strcmp("R-AXIS4   ", magic) == 0)  && (filestat.st_size == 18006000)) {
		printf("Recognized R-AXIS4 100micron...\n");
		/* Copy the frame information */
		*frame_info = staticframeinfo[RAXIS4_100];

		ret_value = 0;
	}
	/* ADSC smv format */
	else if( rdsmv (filename, &smv_header, &lhead, od_pointer, &naxis, axis, &type) == 0) {
		gethd("HEADER_BYTES", smv_value, smv_header);
		lhead = atoi(smv_value);

		fprintf(stderr,"Recognized SMV format. Header = %d bytes. %d x %d pixels\n",
		    lhead,axis[0],axis[1]);

		if (axis[0] == 1152) {
			staticframeinfo[QUANTUM1].header = lhead;
			*frame_info = staticframeinfo[QUANTUM1];
		}
		else {
			staticframeinfo[QUANTUM4].header = lhead;
			*frame_info = staticframeinfo[QUANTUM4];
		}
		close(in_file);
		return(0);
	}
	/* MAR format */
	else if( rdmar (filename, &smv_header, &lhead, od_pointer, &naxis, axis, &type) == 0) {

		fprintf(stderr,"Recognized MAR format. Header = %d bytes. %d x %d pixels\n",
		    lhead,axis[0],axis[1]);

		if (axis[0] == 1200)
			*frame_info = staticframeinfo[MAR180];
		else if(axis[0] == 2000)
			*frame_info = staticframeinfo[MAR300];
		else if(axis[0] == 3000)
			*frame_info = staticframeinfo[MAR345];
		close(in_file);
		return(0);
	}
	/* No format recognized... */
	else {
		printf("Could not recognize format!!\n");
		close(in_file);
		return(2);
	}

	/* Assign memory for the optical densities */
	*od_pointer = realloc(*od_pointer, frame_info->lenx * frame_info->leny * sizeof(od_type));
	if(*od_pointer == NULL) {
		printf("Unable to allocate memory; %d bytes\n",frame_info->lenx * frame_info->leny * sizeof(od_type));
		close(in_file);
		return(3);
	}
	/* Seek past the header  */
	if(frame_info->header > 0)
		lseek(in_file, frame_info->header, SEEK_SET);

	/* read the optical densities */
	read(in_file, *od_pointer, frame_info->lenx * frame_info->leny * sizeof(od_type));

	/* Eventually some post processing */

	/* Apply the correction for the FUJI data */
	if(frame_info->type == FUJI) {
		/* generate the correction table */
		value = dps_malloc("read_frame", "value", 1024*sizeof(od_type));
		for(i=0;i<1024;i++) {
			value[i]= (short) pow(10.0, (double)i/256.0);
		}
		for(i=0;i<frame_info->lenx * frame_info->leny;i++){
			if((*od_pointer)[i] < 1024) {
				(*od_pointer)[i] = value[(*od_pointer)[i]];
			}
			else {
				(*od_pointer)[i] = value[1023];
			}
		}
		free(value);
	}
	min=frame_info->maxod;
	max=0.0;
	for(i=0;i<frame_info->lenx * frame_info->leny;i++){
		if (((*od_pointer)[i] > max ) && ((*od_pointer)[i] < frame_info->overload)) max = (*od_pointer)[i];
		if ((*od_pointer)[i] < min) min = (*od_pointer)[i];
	}

	frame_info->maxod = max;
	frame_info->minod = min;

	close(in_file);
#ifdef READ_FRAME_DEBUG
	printf("READ_FRAME left return value...\n", ret_value);
#endif
	return(ret_value);
}

