/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/


#ifndef BRAVAIS_GOOF_H

#undef BRAVAIS_GOOF_DEBUG

typedef struct {
	float A,B,C,D,E,F;
	} Tensor;


float get_goof(int, struct matrix_3x3 *);  

#define BRAVAIS_GOOF_H
#endif
