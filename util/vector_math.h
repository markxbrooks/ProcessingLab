/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* vector_math.h:  Headerfile for vector_math  */

#ifndef VECTOR_MATH_H

struct dvector {
	double x;
	double y;
	double z;
};

struct vector {
	float x;
	float y;
	float z;
};

struct matrix_3x3 {
	struct vector r[3];
};
 
void cross_prodf(struct vector, struct vector, struct vector *);
float scalar_prodf(struct vector, struct vector);
double scalar_prod(struct dvector, struct dvector);
float vec_lenf(struct vector);
double vec_len(struct dvector);
void cp_vectorf(struct vector, struct vector *);
float vec_volumef(struct vector, struct vector, struct vector);
void invert_3x3(struct matrix_3x3, struct matrix_3x3 *);
void invert_(float *, float *);
void mat_mul_3x3(struct matrix_3x3, struct matrix_3x3, struct matrix_3x3 *);
void mat_add_3x3(struct matrix_3x3, struct matrix_3x3, struct matrix_3x3 *);
void transpose_3x3(struct matrix_3x3, struct matrix_3x3 *);
void print_3x3(struct matrix_3x3);
void copy_3x3(struct matrix_3x3, struct matrix_3x3 *);
void clear_3x3(struct matrix_3x3 *);
void mat_vec_mul_3x3(struct matrix_3x3  m, struct vector v, struct vector *res);
void print_vector_list(int, struct vector *);

#define VECTOR_MATH_H
#endif
