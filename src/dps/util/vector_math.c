/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

#include "vector_math.h"
#include <math.h>
#include <stdio.h>

/* vector_math.c:  Summary of routines for working with vectors and
		   matrices  						*/

/* cross_prodf(): Calcualting cross product of two vectors. Input in    *
 * 		 a and b, output in return				*/

void cross_prodf(struct vector a, struct vector b, struct vector *result)
{
	result->x = (a.y * b.z - b.y * a.z);
	result->y = (a.z * b.x - b.z * a.x);
	result->z = (a.x * b.y - b.x * a.y);
}

/* scalar_prodf(): Calculating scalar product of two vectors. Input in  	*
		  x1,y1,z1 and x2,y2,z2; result in return value		*/

float scalar_prodf(struct vector a, struct vector b)
{
	float prod;

	prod = (a.x) * (b.x) + (a.y) * (b.y) + (a.z) * (b.z);
	return(prod);
}

/* scalar_prod(): Calculating scalar product of two vectors. Input in  	*
		  x1,y1,z1 and x2,y2,z2; result in return value		*
 *		  DOUBLE version					*/

double scalar_prod(struct dvector a, struct dvector b)
{
	double prod;

	prod = (a.x) * (b.x) + (a.y) * (b.y) + (a.z) * (b.z);
	return(prod);
}

/* vec_lenf(): Calculating Calculating the length of a vector on input of*
 *	      a vector a						*/

float vec_lenf(struct vector a)
{
	double tmp;

	tmp = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
	return((float)tmp);
}

/* vec_len(): Calculating Calculating the length of a vector on input of*
 *	      a vector a, DOUBLE version...				*/

double vec_len(struct dvector a)
{
	double tmp;

	tmp = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
	return((float)tmp);
}

/* cp_vectorf(): Copy a given vector a to a vector b 			*/

void cp_vectorf(struct vector a, struct vector *b)

{
	b->x = a.x;
	b->y = a.y;
	b->z = a.z;
}

/* cp_tvectorf(): Copy a given tvector a to a tvector b 			*/

void cp_tvectorf(struct tvector a, struct tvector *b)

{
	b->vec.x = a.vec.x;
	b->vec.y = a.vec.y;
	b->vec.z = a.vec.z;
	b->len = a.len;
	b->fit = a.fit;
}

/* vec_volumef(): Get the volume of a,b and c. 			*/

float vec_volumef(struct vector a, struct vector b, struct vector c)

{
	float volume;   /* obvious :-) */

	volume=fabs(a.x * b.y * c.z +
	    b.x * c.y * a.z +
	    c.x * a.y * b.z -
	    c.x * b.y * a.z -
	    a.x * c.y * b.z -
	    b.x * a.y * c.z);
	return(volume);
}


/* invert_3x3(): Invert 3x3 matrix, using Bob's FORTRAN routine. This	*
 * 		routine should be replaced later by a C-routine		*/

void invert_3x3(struct matrix_3x3  m, struct matrix_3x3 *m_inv)

{
	float m_a[3][3], m_a_inv[3][3]; /* Array as argument for the fortran routine */

	m_a[0][0]=m.r[0].x;
	m_a[1][0]=m.r[1].x;
	m_a[2][0]=m.r[2].x;
	m_a[0][1]=m.r[0].y;
	m_a[1][1]=m.r[1].y;
	m_a[2][1]=m.r[2].y;
	m_a[0][2]=m.r[0].z;
	m_a[1][2]=m.r[1].z;
	m_a[2][2]=m.r[2].z;
	(void)invert_(&m_a[0][0],&m_a_inv[0][0]);
	m_inv->r[0].x=m_a_inv[0][0];
	m_inv->r[1].x=m_a_inv[1][0];
	m_inv->r[2].x=m_a_inv[2][0];
	m_inv->r[0].y=m_a_inv[0][1];
	m_inv->r[1].y=m_a_inv[1][1];
	m_inv->r[2].y=m_a_inv[2][1];
	m_inv->r[0].z=m_a_inv[0][2];
	m_inv->r[1].z=m_a_inv[1][2];
	m_inv->r[2].z=m_a_inv[2][2];
}

/* mat_mul_3x3(): Muliply two 3x3 matrices. Do this by transpose the	*
 *		second matrix and then use the scalar product. 		*/

void mat_mul_3x3(struct matrix_3x3  n, struct matrix_3x3 m, struct matrix_3x3 *out)

{
	struct matrix_3x3 n_t; /* Transposed n */

	transpose_3x3(n, &n_t);

	out->r[0].x=scalar_prodf(m.r[0],n_t.r[0]);
	out->r[0].y=scalar_prodf(m.r[0],n_t.r[1]);
	out->r[0].z=scalar_prodf(m.r[0],n_t.r[2]);
	out->r[1].x=scalar_prodf(m.r[1],n_t.r[0]);
	out->r[1].y=scalar_prodf(m.r[1],n_t.r[1]);
	out->r[1].z=scalar_prodf(m.r[1],n_t.r[2]);
	out->r[2].x=scalar_prodf(m.r[2],n_t.r[0]);
	out->r[2].y=scalar_prodf(m.r[2],n_t.r[1]);
	out->r[2].z=scalar_prodf(m.r[2],n_t.r[2]);
}

void mat_add_3x3(struct matrix_3x3  n, struct matrix_3x3 m, struct matrix_3x3 *out)

{
	out->r[0].x=m.r[0].x+n.r[0].x;
	out->r[0].y=m.r[0].y+n.r[0].y;
	out->r[0].z=m.r[0].z+n.r[0].z;
	out->r[1].x=m.r[1].x+n.r[1].x;
	out->r[1].y=m.r[1].y+n.r[1].y;
	out->r[1].z=m.r[1].z+n.r[1].z;
	out->r[2].x=m.r[2].x+n.r[2].x;
	out->r[2].y=m.r[2].y+n.r[2].y;
	out->r[2].z=m.r[2].z+n.r[2].z;
}

/*transpose_3x3(): Generate the transpose matrix			*/

void transpose_3x3(struct matrix_3x3  m, struct matrix_3x3 *m_t)

{
	m_t->r[0].x=m.r[0].x;
	m_t->r[0].y=m.r[1].x;
	m_t->r[0].z=m.r[2].x;
	m_t->r[1].x=m.r[0].y;
	m_t->r[1].y=m.r[1].y;
	m_t->r[1].z=m.r[2].y;
	m_t->r[2].x=m.r[0].z;
	m_t->r[2].y=m.r[1].z;
	m_t->r[2].z=m.r[2].z;
}

/* print_3x3(): print a 3x3 matrix to spdout			*/

void print_3x3(struct matrix_3x3  a)

{
	printf(" % 10.6f % 10.6f % 10.6f\n",a.r[0].x,a.r[1].x,a.r[2].x);
	printf(" % 10.6f % 10.6f % 10.6f\n",a.r[0].y,a.r[1].y,a.r[2].y);
	printf(" % 10.6f % 10.6f % 10.6f\n",a.r[0].z,a.r[1].z,a.r[2].z);
}
/* fprint_3x3(): print a 3x3 matrix to spdout			*/

void fprint_3x3(char *filename, struct matrix_3x3  a)

{
	FILE *output;
	output = fopen(filename, "a");

	fprintf(output,"\n % 10.6f % 10.6f % 10.6f\n",a.r[0].x,a.r[1].x,a.r[2].x);
	fprintf(output," % 10.6f % 10.6f % 10.6f\n",a.r[0].y,a.r[1].y,a.r[2].y);
	fprintf(output," % 10.6f % 10.6f % 10.6f\n",a.r[0].z,a.r[1].z,a.r[2].z);
}

/*copy_3x3(): Generate a copy of a matrix			*/

void copy_3x3(struct matrix_3x3  in, struct matrix_3x3 *out)

{
	out->r[0].x=in.r[0].x;
	out->r[0].y=in.r[0].y;
	out->r[0].z=in.r[0].z;
	out->r[1].x=in.r[1].x;
	out->r[1].y=in.r[1].y;
	out->r[1].z=in.r[1].z;
	out->r[2].x=in.r[2].x;
	out->r[2].y=in.r[2].y;
	out->r[2].z=in.r[2].z;
}

/*clear_3x3(): Generate a copy of a matrix			*/

void clear_3x3(struct matrix_3x3  *om)

{
	om->r[0].x=0.0;
	om->r[0].y=0.0;
	om->r[0].z=0.0;
	om->r[1].x=0.0;
	om->r[1].y=0.0;
	om->r[1].z=0.0;
	om->r[2].x=0.0;
	om->r[2].y=0.0;
	om->r[2].z=0.0;
}

/*mat_vec_mul_3x3(): Multiply vector by matrix resulting in vector	*/

void mat_vec_mul_3x3(struct matrix_3x3  m, struct vector v, struct vector *res)

{
	res->x=v.x*m.r[0].x+v.y*m.r[1].x+v.z*m.r[2].x;
	res->y=v.x*m.r[0].y+v.y*m.r[1].y+v.z*m.r[2].y;
	res->z=v.x*m.r[0].z+v.y*m.r[1].z+v.z*m.r[2].z;
}

/*add_vecf(): Add 2 vectors	*/

void add_vecf(struct vector  a, struct vector b, struct vector *res)

{
	res->x=a.x+b.x;
	res->y=a.y+b.y;
	res->z=a.z+b.z;
}

/* print_vector_list(): Prints a list of vectors columnwise.		*/

void print_vector_list(int n, struct vector v[])

{
	int i;		/* Counter */

	printf("     X          Y          Z\n");
	for(i=0;i<n;i++) {
		printf("% 10.8f % 10.8f % 10.8f\n",v[i].x,v[i].y,v[i].z);
	}
}

void abc_to_cell(float *abc, struct matrix_3x3 *m) {
	m->r[0].x = abc[0]; m->r[0].y = 0; m->r[0].z = 0;
	m->r[1].x = abc[1] * cos(M_PI / 180.0 * abc[5]);
	m->r[1].y = abc[1] * sin(M_PI / 180.0 * abc[5]);
	m->r[1].z = 0;
	m->r[2].x = abc[2] * cos(M_PI / 180.0 * abc[4]);
	m->r[2].y = abc[2] * (cos(M_PI / 180.0 * abc[3]) - 
			      cos(M_PI / 180.0 * abc[5]) *
			      cos(M_PI / 180.0 * abc[4])) / sin(M_PI / 180.0 * abc[5]);
	m->r[2].z = sqrt(abc[2] * abc[2] * sin(M_PI / 180.0 * abc[4]) *
			 sin(M_PI / 180.0 * abc[4]) - m->r[2].y * m->r[2].y);
}

float vec_angle(struct vector a, struct vector b) {
	float len_a = vec_lenf(a), len_b = vec_lenf(b);
	if (len_a * len_b < 0.000001) return 0.0;
	return acos(scalar_prodf(a, b) / len_a / len_b);
}

void cell_to_abc(struct matrix_3x3 *m, float *abc) {
	abc[0] = vec_lenf(m->r[0]);
	abc[1] = vec_lenf(m->r[1]);
	abc[2] = vec_lenf(m->r[2]);
	abc[3] = vec_angle(m->r[1], m->r[2]) / M_PI * 180.0;
	abc[4] = vec_angle(m->r[0], m->r[2]) / M_PI * 180.0;
	abc[5] = vec_angle(m->r[0], m->r[1]) / M_PI * 180.0;
}
