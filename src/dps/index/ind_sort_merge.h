/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* ind_sort_merge.h:  Header File for ind_sort_merge.c  */

#ifndef IND_SORT_MERGE

#include "index_def.h"

#undef MERGE_TVECTORS_SAME_N_DEBUG
#undef CUT_TVECTOR_LIST

#define MERGE_DIFF_ABS	1.0	/* Maximum ABSOLUTE difference between
				   t-vectors for merging		*/

#define MERGE_DIFF_REL	0.1	/* Maximum RELATIVE difference between
				   t-vectors for merging		*/

#define CUT_CRIT	0.05	/* Criterium, how much the fit value is 
				   allowed to drop before the tvector
				   list is cut				*/

void sort_tvectors_on_n(int, struct tvector *);
void sort_tvectors_on_l(int, struct tvector *);
void sort_tvectors_on_n_and_l(int, struct tvector *);
int merge_tvectors_same_n(int, struct tvector *);
int cut_tvector_list(int, float, struct tvector *);
int same_tvector (struct tvector , struct tvector , float, int *);
#define IND_SORT_MERGE

#endif
