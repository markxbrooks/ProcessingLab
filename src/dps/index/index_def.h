/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* index.h:  Headerfile for index  */

#ifndef INDEX_DEF_H

#include "util.h"

#undef INDEX_DEBUG
#undef READINPUT_DEBUG
#undef READINFO_DEBUG
#undef READ_RS_DEBUG
#undef READXYZ_DEBUG

#if defined __linux__ || defined (__APPLE__) || defined (__sun)
#include "linux_defs.h"
#endif /* linux or apple or sun */

#if defined _WIN32 || defined (_MSC_VER) || defined (WIN32) || defined (__WIN32)
#include "windows_defs.h"
#endif /* linux or apple or sun */

#define N_TVECTOR 150		/* maximum number of generated t-vectors */
                                /* Note that actual number used is passed from MOSFLM */

#define MAX_OM_COUNT 30		/* Maximum number of orientation matrices that 
				   will be generated from t_vectors and checked */

typedef short Index_Mode;	/* Type for the mode of indexing */
#define INDEX_AUTOMATIC 0


/* Type for defining the Mode, the indexing routine should run. */

int do_index(int, int, int, int, float, int, float, int, int, int, int, int, float, float,
             float, float, float, float, int *, int, int, int, int,
             int, struct vector *, float, float, float*, Index_Mode, struct matrix_3x3 *);  

#define INDEX_DEF_H
#endif
