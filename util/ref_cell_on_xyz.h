#ifndef REF_CELL_ON_XYZ_H

#define CELL_REF_MAX_DIFF_HKL 0.2
#define CELL_N_PAR 6
#define OM_N_PAR 9

void ref_cell_on_xyz(int number_xyz, struct vector *xyz, struct matrix_3x3 *om);
void ref_om_on_xyz(int number_xyz, struct vector *xyz, struct matrix_3x3 *om);

#define REF_CELL_ON_XYZ_H
#endif
