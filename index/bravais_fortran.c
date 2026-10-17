#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../util/util.h"
#include "bravais.h"

/* bravais_fortran_(): This routine converts array from FORTRAN program
                       to structure for call to bravais
		       31.03.2010 now needs and returns orientation matrices
		       in Mosflm orientation

   bravais() is the Cambridge 5-argument form (see index/bravais.c).
   Fortran still receives *error; set to 0 after a successful call.
*/

void bravais_fortran_(float *a, float *wave, struct matrix_3x3 cl_om[45], float crit[45], char sys_id[45], int *error)  
{
    int i, output_type;
    struct matrix_3x3 dpstomosflm;
    struct matrix_3x3 om;           /* Generated orientation matrices */
    om.r[0].z = a[0]/(*wave);
    om.r[0].x = a[1]/(*wave);
    om.r[0].y = a[2]/(*wave);
    om.r[1].z = a[3]/(*wave);
    om.r[1].x = a[4]/(*wave);
    om.r[1].y = a[5]/(*wave);
    om.r[2].z = a[6]/(*wave);
    om.r[2].x = a[7]/(*wave);
    om.r[2].y = a[8]/(*wave);
    output_type = 1;
    bravais(&om, cl_om, crit, sys_id, output_type);
    if (error != NULL) {
        *error = 0;
    }
    for (i = 0; i < 45; i++ ){
      dpstomosflm.r[0].x = *wave * cl_om[i].r[0].z;
      dpstomosflm.r[0].y = *wave * cl_om[i].r[0].x;
      dpstomosflm.r[0].z = *wave * cl_om[i].r[0].y;
      dpstomosflm.r[1].x = *wave * cl_om[i].r[1].z;
      dpstomosflm.r[1].y = *wave * cl_om[i].r[1].x;
      dpstomosflm.r[1].z = *wave * cl_om[i].r[1].y;
      dpstomosflm.r[2].x = *wave * cl_om[i].r[2].z;
      dpstomosflm.r[2].y = *wave * cl_om[i].r[2].x;
      dpstomosflm.r[2].z = *wave * cl_om[i].r[2].y;
      copy_3x3(dpstomosflm, &cl_om[i]);
    }
    return;
}


