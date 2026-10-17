/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* ind_choose.h:  Headerfile for ind_choose.f */

#undef CHOOSE_FIRST_NON_LINEAR_DEBUG
#undef CHOOSE_BEST_OM_DEBUG


/* The next two definitions are for the routine choose_best_om(). Usage
   see there.								*/

#define VOLUME_FACTOR   1.6  /* was 1.6 */
#define MAX_DIFF_HKL	0.30 /* was 0.2 */

/* Definition for choose_non_lin(). 					*/

#define LIN_DEP_FACTOR 0.5	/* was 0.2 Factor between volume and product of
				   the vector length for the decision of
				   linear independence			*/

int choose_non_lin(int, float, float, float, float, int, struct vector xyz[], int, struct tvector t[], int, float *, struct matrix_3x3 *);
int choose_non_lin_bf(int, float, float, float, float, int, struct vector xyz[], int, struct tvector t[], int, float *, struct matrix_3x3 *);
int choose_best_om(int, float, int, struct matrix_3x3 *, int, struct vector *, int);

