#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../util/util.h"

/* write_ub_matrix(): This routine uses DPS routines to
                    to print some bits and bobs from UB matrix.
*/

void write_ub_matrix_(float *a, float u[9], float mis_set[3])  

{
  struct matrix_3x3 om,bmatrix,umatrix;
  struct setangles_struct setangles;
  struct cell_struct cell;
  float rtod = 180.0/3.141592654;

  om.r[0].x = a[0];
  om.r[0].y = a[1];
  om.r[0].z = a[2];
  om.r[1].x = a[3];
  om.r[1].y = a[4];
  om.r[1].z = a[5];
  om.r[2].x = a[6];
  om.r[2].y = a[7];
  om.r[2].z = a[8];
  
  dircell_from_om(om,&cell);
  
  B_and_U_from_om(om,&bmatrix,&umatrix);
  
  u[0] = umatrix.r[0].x;
  u[1] = umatrix.r[0].y;
  u[2] = umatrix.r[0].z;
  u[3] = umatrix.r[1].x;
  u[4] = umatrix.r[1].y;
  u[5] = umatrix.r[1].z;
  u[6] = umatrix.r[2].x;
  u[7] = umatrix.r[2].y;
  u[8] = umatrix.r[2].z;

  setangles_from_reccell(cell,om,&setangles);

  mis_set[0] = setangles.x;
  mis_set[1] = setangles.y;
  mis_set[2] = setangles.z;
}


