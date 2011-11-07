/**********************************************************************
 *
 * npem
 * version 0.50
 * 22 June 2003
 *
 * Karl W Broman
 *
 * npem_matrix.c : ludcmp, lubksb, invert_matrix
 *
 * 4/26/96 (revised 5/3/96, 8/18/2000, 6/22/2003)
 *
 * These are a set of routines to do a matrix inversion.
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "npem.h"
#include <R.h>
#include <R_ext/PrtUtil.h>

/**********************************************************************
 *
 * ludcmp
 *
 *      This is a routine from "Numerical Recipes in C" to perform
 *  an LU decomposition.  (See pages 46-47.)
 *
 *  INPUT:
 *
 *         a = (n x n) matrix on which to perform the decomposition.
 *             The matrix is destroyed, being replaced by the LU decom.
 *
 *         n = size of matrix (num rows = num cols = n)
 *
 *         index = output vector (of length n) recording the row
 *                 permutation effected by the partial pivoting.
 *
 *         d = output as +/- 1 depending on whether the number of row
 *             interchanges is even or odd, respectively.
 *
 *         work = vector of doubles, of length n, to use as workspace.
 *
 *         errorflag = output as 1 if error, 0 otherwise
 *
 **********************************************************************/

void ludcmp(double *a, int *n, int *indx, double *d, double *work,
            int *errorflag)
{
  int i, imax, j, k;
  double big, dum, sum, temp;

  *d = 1.0;
  *errorflag = 0;
  for(i=0; i< *n; i++) {
    big = 0.0;
    for(j=0; j< *n; j++)   /* loop over rows to get scaling information */
      if((temp=fabs(a[i+ *n * j])) > big) big=temp;
    if(big == 0.0) {
      Rprintf("Singular matrix in routine ludcmp.\n");
      *errorflag = 1;
      return;
    }
    work[i] = 1.0/big;   /* save the scaling */
  }

  for(j=0; j< *n; j++) { /* loop of columns of Crout's method */
    for(i=0; i<j; i++) {
      sum = a[i +  *n * j];
      for(k=0; k<i; k++) sum -= (a[i+ *n * k] * a[k+ *n * j]);
      a[i+ *n * j] = sum;
    }
    big = 0.0;   /* initialize for the search for largest pivot element */
    for(i=j; i< *n; i++) {
      sum=a[i+ *n *j];
      for(k=0; k<j; k++) sum -= (a[i+ *n *k] * a[k+ *n * j]);
      a[i+ *n * j] = sum;
      if((dum=work[i]*fabs(sum)) >= big) {
        big = dum;   /* is the figure of merit for the pivot better */
        imax = i;    /*      than the best so far? */
      }
    }
    if(j != imax) {  /* do we need to interchange rows? */
      for(k=0; k< *n; k++) {      /* yes */
        dum = a[imax+ *n * k];
        a[imax+ *n * k] = a[j + *n * k];
        a[j + *n * k] = dum;
      }
      *d = -(*d);          /* change parity of d */
      work[imax] = work[j];   /* and interchange scale factor */
    }
    indx[j] = imax;
    if(a[j+ *n *j] == 0.0) a[j+ *n *j] = 1.0e-20;

    if(j != (*n-1)) {  /* now, finally, divide by the pivot element */
      dum = 1.0/a[j+ *n*j];
      for(i=j+1; i< *n; i++) a[i+ *n * j] *= dum;
    }
  }
}




/**********************************************************************
 *
 * lubksb:
 *      This is a routine from "Numerical Recipes in C" to perform
 *  forward and backward substitution using the LU decomposition of
 *  a matrix.  Thus, in conjunction with "ludcmp," it solves the
 *  matrix equation AX = b.  (See page 47.)
 *
 *  INPUT:
 *
 *         a = (n x n) matrix, output from "ludcmp"
 *
 *         n = size of matrix (num rows = num cols = n)
 *
 *         index = vector of length n, output from "ludcmp"
 *
 *         b = vector of length n (right hand side of above equation)
 *
 **********************************************************************/

void lubksb(double *a, int *n, int *indx, double *b)
{
  int i, ii= -1, ip, j;
  double sum;

  for(i=0; i< *n; i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii >= 0) for(j=ii; j<=i-1; j++) sum -= a[i+ *n * j] *b[j];
    else if(sum) ii=i;
    b[i]=sum;
  }
  for(i= *n - 1; i>=0; i--) {
    sum = b[i];
    for(j=i+1; j < *n; j++) sum -= a[i+ *n*j]*b[j];
    b[i] = sum/a[i+ *n*i];
  }
}



/**********************************************************************
 *
 * invert_matrix
 *
 *      This routine inverts a matrix, using the above routines from
 *  "Numerical Recipes in C" (see page 48).
 *
 * INPUT:
 *
 *       a = n x n matrix to be inverted.  The contents are destroyed.
 *
 *       y = n x n matrix in which inverse is to be placed.
 *
 *       n = size of matrix
 *
 *       intwork = vector of length n, to be used as workspace
 *
 *       doublework = vector of length n, to be used as workspace
 *
 *       errorflag = output as 1 if error, 0 otherwise
 *
 **********************************************************************/

void invert_matrix(double *a, double *y, int *n, int *intwork,
                   double *doublework, int *errorflag)
{
  double d;
  int i, j;

  ludcmp(a, n, intwork, &d, doublework, errorflag);
  if(*errorflag) return;

  for(j=0; j< *n; j++) {
    for(i=0; i < *n; i++) doublework[i]=0.0;
    doublework[j] = 1.0;
    lubksb(a, n, intwork, doublework);
    for(i=0; i< *n; i++) y[i+ *n * j] = doublework[i];
  }
}
