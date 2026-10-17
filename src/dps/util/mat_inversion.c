/*=====================================================================*
 *                                                                     *
 *                         Data Processing Suit                        *
 *                                                                     *
 *                            Utility-Library                          *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                         File: mat_inverion.c                        *
 *                                                                     *
 *=====================================================================*/

/* Thsi is still the matrix inversion routine from 'Numerical Recipies'.
   Should be replaced as soon as possible......			       */

#include <math.h> 
#include "nrutil.h" 
#include "mat_inversion.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}


/* gaussj() matrix inversion after Gauss-Jordan, giving a one-dimensional
   array holding the matrix to be inverted; b holds the right-hand
   sides....								*/

void gaussj(double *a, int n, float *b, int m)

/* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. 
a[1..n][1..n] is the input matrix. b[1..n][1..m] is input containing the m 
right-hand side vectors. On output, a is replaced by its matrix inverse, and 
b is replaced by the corresponding set of solution vectors. */

{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	/* The integer arrays ipiv, indxr, andindxc are used for bookkeeping 
	   on the pivoting. */
	indxc=ivector(1,n-1);
	indxr=ivector(1,n-1);
	ipiv=ivector(1,n-1);
	for (j=0;j<n;j++) ipiv[j]=0;

	/* This is the main loop over the columns to be reduced. */

	for (i=0;i<n;i++) {
		big=0.0;
		/* This is the outer loop of the search for a pivot element. */

		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j*n+k]) >= big) {
							big=fabs(a[j*n+k]);
							irow=j; 
							icol=k;
						}
					}
					else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);

		/* We now have the pivot element, so we interchange rows, if needed, 
to put the pivot element on the diagonal. The columns are not physically 
interchanged, only relabeled: indxc[i], the column of the ith pivot element, 
is the ith column that is reduced, while indxr[i] is the row in which that 
pivot element was originally located. If indxr[i] 6 = indxc[i] there is an 
implied column interchange. With this form of bookkeeping, the solution b's 
will end up in the correct order, and the inverse matrix will be scrambled by 
columns. */
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow*n+l],a[icol*n+l])
				for (l=0;l<m;l++) SWAP(b[irow*m+l],b[icol*m+l])
		}
		/* We are now ready to divide the pivot row by the pivot element, located at 
	   irow and icol. */

		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol*n+icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol*n+icol];
		a[icol*n+icol]=1.0;
		for (l=0;l<n;l++) {
			a[icol*n+l] *= pivinv;
		}
		for (l=0;l<m;l++) {
			b[icol*m+l] *= pivinv;
		}
		/* Next, we reduce the rows... */

		for (ll=0;ll<n;ll++)
			/* ...except for the pivot one, of course. */

			if (ll != icol) {
				dum=a[ll*n+icol];
				a[ll*n+icol]=0.0;
				for (l=0;l<n;l++) {
					a[ll*n+l] -= a[icol*n+l]*dum;
				}
				for (l=0;l<m;l++) {
					b[ll*m+l] -= b[icol*m+l]*dum;
				}
			}
	}

	/* This is the end of the main loop over columns of the reduction. 
It only remains to unscram- ble the solution in view of the column 
interchanges. We do this by interchanging pairs of columns in the reverse 
order that the permutation was built up. */

	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++) SWAP(a[k*n+indxr[l]],a[k*n+indxc[l]]);
	}

	/* And we are done. */

	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

/* gaussj() matrix inversion after Gauss-Jordan, giving a two-dimensional
   array in C holding the matrix to be inverted; b holds the right-hand
   sides....								*/

void gaussj_mat(double **a, int n, float *b, int m)

/* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. 
a[1..n][1..n] is the input matrix. b[1..n][1..m] is input containing the m 
right-hand side vectors. On output, a is replaced by its matrix inverse, and 
b is replaced by the corresponding set of solution vectors. */

{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	/* The integer arrays ipiv, indxr, andindxc are used for bookkeeping 
	   on the pivoting. */
	indxc=ivector(1,n-1);
	indxr=ivector(1,n-1);
	ipiv=ivector(1,n-1);
	for (j=0;j<n;j++) ipiv[j]=0;

	/* This is the main loop over the columns to be reduced. */

	for (i=0;i<n;i++) {
		big=0.0;
		/* This is the outer loop of the search for a pivot element. */

		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j; 
							icol=k;
						}
					}
					else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);

		/* We now have the pivot element, so we interchange rows, if needed, 
to put the pivot element on the diagonal. The columns are not physically 
interchanged, only relabeled: indxc[i], the column of the ith pivot element, 
is the ith column that is reduced, while indxr[i] is the row in which that 
pivot element was originally located. If indxr[i] 6 = indxc[i] there is an 
implied column interchange. With this form of bookkeeping, the solution b's 
will end up in the correct order, and the inverse matrix will be scrambled by 
columns. */
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
				for (l=0;l<m;l++) SWAP(b[irow*m+l],b[icol*m+l])
		}
		/* We are now ready to divide the pivot row by the pivot element, located at 
	   irow and icol. */

		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) {
			a[icol][l] *= pivinv;
		}
		for (l=0;l<m;l++) {
			b[icol*m+l] *= pivinv;
		}
		/* Next, we reduce the rows... */

		for (ll=0;ll<n;ll++)
			/* ...except for the pivot one, of course. */

			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) {
					a[ll][l] -= a[icol][l]*dum;
				}
				for (l=0;l<m;l++) {
					b[ll*m+l] -= b[icol*m+l]*dum;
				}
			}
	}

	/* This is the end of the main loop over columns of the reduction. 
It only remains to unscram- ble the solution in view of the column 
interchanges. We do this by interchanging pairs of columns in the reverse 
order that the permutation was built up. */

	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++) SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}

	/* And we are done. */

	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

/* invert_2x2(): Inverting a 2 by 2 matrix as an array by using the determinant
    ant the following formula...			
			[a22 -a12]
	A-1 = 1/(det A) [	 ]
			[-a21 a11]
									*/
void invert_2x2(float in[2][2], float out[2][2])

{
	float det; /* Determinante */

	det=in[0][0]*in[1][1]-in[0][1]*in[1][0];
	out[0][0]=in[1][1]/det;
	out[1][0]=-in[1][0]/det;
	out[0][1]=-in[0][1]/det;
	out[1][1]=in[0][0]/det;
}
