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
 *                           File: allocate.c                          *
 *                                                                     *
 *=====================================================================*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#if !defined vms && !defined(__APPLE__) && ! ( defined (__FreeBSD__) && defined (__STDC__) )
#include <malloc.h>
#endif

/* dps_malloc(): Memory allocation using malloc with error handling.
 *               The procedure exits if a momory allocation error 
 *               occures 						*/

void *dps_malloc( char *routine, char *variable, size_t size)

{
	void *pointer;

	pointer = malloc(size);
	if (pointer == NULL) {
		fprintf(stderr,"\n Fatal memory allocation error in routine %s for variable %s of size %d bytes\n",routine,
		    variable,(int)size);
		exit(1);
	}
	return(pointer);
}
