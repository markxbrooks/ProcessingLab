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
 *                              dps_bravais                            *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                        File: dps_bravais_read.c                     *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "util.h"
#include "bravais.h"
#include "dps_bravais_read.h"

/* read_control(): Read the Control input file and analyse it. Return the
    		  control structure and an exit value for successfull
		  reading.						*/

Boolean read_control_do_read(FILE *in_file,  bravais_control *control)

{
	int i;
	char line[256];
	char command[100];
	char *rest;
	char name[100];
	char tmp_mode[20];
	float start, end;
	int line_count = 0;
	Boolean file_ok=TRUE;
	double pi = M_PI;
	Boolean command_identified;

#define MAX_PEAKFILES 2000

#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL DO_READ entered...\n");
#endif

	while (fgets(line, sizeof(line), in_file) != 0) {
		command_identified = FALSE;
		line_count++;
		if(sscanf(line, "%s", command) <= 0) {
			command[0]='\0';
			command_identified = TRUE;
		}
		if(command[0] == '#') {
			command[0]='\0';
			command_identified = TRUE;
		}
		rest = &(line[strlen(command)]);
		upper(command);
#ifdef READ_CONTROL_DEBUG
		printf("line %d: %s\n", line_count, command);
#endif
		if (strcmp(command, "OMFILE") == 0) {
			command_identified = TRUE;
			switch(sscanf(rest, "%s", command)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the control input file: No data after ORIENTATIONMATRIX keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 1: 
				control->om_filename = makestring(command);
				break;
			default:
				printf("Error in line %d in the control input file: To many parameters after ORIENTATIONMATRIX keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "MODE") == 0) {
			command_identified = TRUE;
			if(sscanf(rest, "%s", command) <= 0) command[0]='\0';
			upper(command);
			if(strcmp(command, "COMPARE") == 0) {
				control->mode=BRAVAIS_COMPARE;
			}
			else control->mode=BRAVAIS_DETERMINE;
		}
		if (command_identified == FALSE) {
			printf("Error in line %d in the control input file: Keyword %s unknown\n", line_count, command);
			file_ok=FALSE;
		}
	}
	if((control->mode == BRAVAIS_COMPARE) && (control->om_filename[0] == '\0')) {
		printf("ORIENTATIONMATRIX needed if MODE = COMPARE!!!\n");
		file_ok=FALSE;
	}
#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL DO_READ left...\n");
#endif
	return(file_ok);
}


Boolean read_control(char *filename, bravais_control *control)

{

	FILE *in_file;
	Boolean status;

#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL entered...\n");
#endif
	in_file = fopen(filename, "r");
	if (in_file == NULL) {
		printf("Can't open %s!! Exit/n", filename);
		return(FALSE);
	};
	status = read_control_do_read(in_file, control);
	fclose(in_file);
#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL left...\n");
#endif
	return(status);

}


