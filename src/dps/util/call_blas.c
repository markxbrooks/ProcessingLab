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
 *                           File: call_blas.c                         *
 *                                                                     *
 *=====================================================================*/

/* call_blas.c:  Summary of routines for calling BLAS routines from C	*/
/* ATTENTION: THE ROUTINES IN THIS FILE ARE NOT USED AT THE MOMENT. IT
              CANNOR BE GUARANTEED THAT THEY ARE WORKING PROPERLY.      */
#include "call_blas.h"

/* mat_mul_matT(): Generate the Product of a matrix and it's transpose	*
 *		using sgemm_ (BLAS3).					*/

void mat_mul_matT(float *a, int m, int n, float *res)
/* Arguments are: a: n x m (n rows, m coloumns) matrix (def: a[n][m]!!! 
		  Let n be the fast running index!!! C -> FORTRAN Probelm!!!)
		  res: ergebnis matrix n x n				*/
{
	char n_ch = 'n';
	char t_ch = 't';
	float one = 1.0;
	float zero = 0.0;

	sgemm_(&t_ch,&n_ch,&n,&n,&m,&one,a,&m,a,&m,&zero,res,&n);
}

/* matT_mul_vec(): Generate the Product of the transpose of a matrix	*
 *		and a vector using sgemm_ (BLAS3).			*/

void matT_mul_vec(float *a, int m, int n, float *x, float *y)
/* Arguments are: a: n x m (n rows, m coloumns) matrix (def: a[n][m]!!! 
		  Let n be the fast running index!!! C -> FORTRAN Problem!!!)
		  x: input vector, y: output vector 			*/
{
	char n_ch = 'n';
	char t_ch = 't';
	float one = 1.0;
	float zero = 0.0;

	sgemv_(&t_ch,&n,&m,&one,a,&m,x,&one,&zero,y,&one);
}
