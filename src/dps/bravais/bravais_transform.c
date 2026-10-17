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
 *                       File: bravais_transform.c                     *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "bravais.h"
#include "bravais_transform.h"

/* MAtrix with the transformation matrices for all lattice characters */

struct matrix_3x3 tm[45] = {
	{ 0,  0,  0,  0,  0,  0,  0,  0,  0},
	{ 1, -1,  1,  1,  1, -1, -1,  1,  1},
	{ 1, -1,  0, -1,  0,  1, -1, -1, -1},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1},
	{ 1, -1,  0, -1,  0,  1, -1, -1, -1},
	{ 1,  0,  1,  1,  1,  0,  0,  1,  1},
	{ 0,  1,  1,  1,  0,  1,  1,  1,  0},
	{ 1,  0,  1,  1,  1,  0,  0,  1,  1},
	{-1, -1,  0, -1,  0, -1,  0, -1, -1},
	{ 1,  0,  0, -1,  1,  0, -1, -1,  3},
	{ 1,  1,  0,  1, -1,  0,  0,  0, -1},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1},
	{ 1,  1,  0, -1,  1,  0,  0,  0,  1},
	{ 1,  1,  0, -1,  1,  0,  0,  0,  1},
	{ 1,  0,  0,  0,  1,  0,  1,  1,  2},
	{-1, -1,  0,  1, -1,  0,  1,  1,  2},
	{ 1, -1,  0, -1, -1,  0, -1,  0, -1},
	{ 0, -1,  1,  1, -1, -1,  1,  0,  0},
	{-1,  0,  0,  0, -1,  1, -1,  1,  1},
	{ 0,  1,  1,  0,  1, -1, -1,  0,  0},
	{ 0,  1,  0,  0,  0,  1,  1,  0,  0},
	{ 0,  1,  0,  0,  0,  1,  1,  0,  0},
	{ 0,  1,  1,  0, -1,  1,  1,  0,  0},
	{ 1,  2,  1,  0, -1,  1,  1,  0,  0},
	{ 0,  1,  1,  0, -1,  1,  1,  0,  0},
	{ 1,  0,  0, -1,  2,  0, -1,  0,  2},
	{-1,  2,  0, -1,  0,  0,  0, -1,  1},
	{-1,  0,  0, -1,  0,  2,  0,  1,  0},
	{ 1,  0,  0,  1, -2,  0,  0,  0, -1},
	{ 0,  1,  0,  0,  1, -2, -1,  0,  0},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1},
	{-1,  0,  0,  0,  0, -1,  0, -1,  0},
	{ 0, -1,  0, -1,  0,  0,  0,  0, -1},
	{ 1,  0,  0, -1,  0, -2,  0,  1,  0},
	{ 1,  0,  2,  1,  0,  0,  0,  1,  0},
	{-1,  0,  0,  1,  2,  0,  0,  0, -1},
	{-1, -2,  0, -1,  0,  0,  0,  0, -1},
	{ 0, -1,  0,  0,  1,  2, -1,  0,  0},
	{ 0, -1, -2,  0, -1,  0, -1,  0,  0},
	{-1,  0,  0,  0, -1,  0,  1,  1,  2},
	{-1,  0,  0, -1, -1, -2,  0, -1,  0},
	{ 1,  0,  0,  0,  1,  0,  0,  0,  1}
};
/* 44 aP */

void transform_om(int lattice, struct matrix_3x3 *om_in, struct matrix_3x3 *om_out)

{
	mat_mul_3x3(*om_in, tm[lattice], om_out);
}


