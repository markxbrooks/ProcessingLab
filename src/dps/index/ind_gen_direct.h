/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* ind_gen_direct.h:  Header File to ind_gen_direct.h  */

#ifndef IND_GEN_DIRECT

#include "util.h"
#include "index_def.h"

#undef GEN_PROJECTIONS_DEBUG
#undef GEN_DIRECT_LATTICE_VECTORS_DEBUG
#define GET_FITTING_DEBUG 0
/* make this = 2 to get fft & projection files output */
#define NEW_GET_FITTING_DEBUG 0
#undef REF_N_VECTOR_DEBUG
/* hrp 24.10.2006 write out FFT first big peaks to a big 2D array */
#undef WRITE_2D_ARRAY
#undef WRITE_REVERSE
#define REVERSE_SIZE 201
struct group_struct {
	int start,end;
	float mean;
	struct vector n_vector;
	float chi_sqr;
};

struct ds_struct {
	int start,end;
	int num;
};

int  gen_direct_lattice_vectors(int, int, int, int, float, int, float, 
             int, int,int,int, int,float,float,int *, int, int,
             int, float, float, struct vector *, int *,struct tvector *);

void rffti_(int *,float *);
void rfftf_(int *,float *,float *);

#define IND_GEN_DIRECT

#endif
