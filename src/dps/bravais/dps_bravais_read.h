/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* dps_index_read.h:  Header File to dps_index_read.c  */

#ifndef DPS_BRAVAIS_READ

#undef READ_CONTROL_DEBUG

#include "util.h"
#include "bravais.h"

typedef struct {
        BravaisMode      mode;   	/* Mode */
        char	 	 *om_filename;	/* Name of OM to be compared */
} bravais_control;


Boolean read_control(char *, bravais_control *);
Boolean read_control_do_read(FILE *, bravais_control *);

#define DPS_BRAVAIS_READ

#endif
