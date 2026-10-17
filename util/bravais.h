/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/


#ifndef BRAVAIS_H

#include "bravais_help.h"
#include "bravais_goof.h"
#include "bravais_permutate.h"
#include "bravais_transform.h"

#define BRAVAIS_DEBUG

typedef short BravaisMode;

#define BRAVAIS_COMPARE 0
#define BRAVAIS_DETERMINE 1

/* Match index/bravais.c and src/dps/bravais/bravais.c (5-arg API). */
void bravais(struct matrix_3x3 *, struct matrix_3x3 *, float crit[45], char sys_id[45][3], int output_type);
#define BRAVAIS_H
#endif
