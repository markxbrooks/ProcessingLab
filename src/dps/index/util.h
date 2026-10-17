#ifndef UTIL_H


#include "util_def.h"
#include "vector_math.h"
#include "cryst_matrix.h"
#include "read_peaks.h"
#include "ind_transform.h"
#include "nrutil.h"
#include "mat_inversion.h"
#include "call_blas.h"
#include "ind_stat.h"
#include "fftn.h"
#include "read_om.h"
#include "allocate.h"
#include "ref_cell_on_xyz.h"
#include "read_frame.h"
char *makestring( char *s );
void upper(char *s);

#define UTIL_H
#define RS_2_XYZ_DEBUG TRUE
#endif
