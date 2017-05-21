/**********************************************************************
 *
 * npem
 * version 0.50
 * 22 June 2003
 *
 * Karl W Broman
 *
 * npem_start.c : npem_start, piksrt
 *
 * 3/13/96 (revised 4/25/96, 8/20/2000, 6/22/2003)
 *
 * See the file "npem.h" for a thorough description.
 *
 **********************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"npem.h"
#include <R.h>
#include <R_ext/PrtUtil.h>

/**********************************************************************
 *
 * npem_start
 *
 * This function finds a crude starting point for the EM algorithm;
 * it is not completely general.  Each plate should contain only one
 * cell density.
 *
 * Input: (all are pointers to either double or int)
 *        The input is for the most part a subset of that for npem_em,
 *        except for "n_sd" and "work"
 *
 *      y           the vector of scintillation counts y{pij}.
 *                  These should be in "lexicographic order,"
 *                  with the wells moving most quickly, then antigen
 *                  groups, and finally plates.  All plates should
 *                  contain at least one well for each block.
 *
 *      n_plates    the number of plates in the NPEM
 *
 *      c           the vector of cell densities c{pij}.
 *                  c{pij} gives the cell density corresponding to
 *                  y{pij}.  (Usually this is either a vector of 1's
 *                  or cells/well / 10^6 cells.)
 *
 *      n_groups    the number of groups on each plate
 *
 *      n           n{pi} = number of wells for group i on plate p
 *                  (a vector of length n_plates x n_groups)
 *
 *        ns and len_plate can be empty on input,
 *      they'll be calculated from n
 *
 *      ns          ns{pi} = partial sums n = (0,n[0],n[0]+n[1],...)
 *                  (same length as n)
 *
 *      len_plate   len_plate{p} = number of wells used on plate p
 *                  (vector of length n_plates) = sum_i n{pi}
 *
 *      cv          common coefficient of variation (= SD/ave), used in
 *                  randomizing the starting point; if cv=0, don't
 *                  randomize
 *
 *      n_sd        number of SDs to use as cut-off
 *
 *      ests        will take initial estimates (lambda's, (a, b, sigma)'s)
 *
 *      work        workspace of doubles, of length
 *                      n_plates (2 n_groups + 2)
 *
 **********************************************************************/

void npem_start(double *y, int *n_plates, double *c, int *n_groups,
               int *n, int *ns, int *len_plate, double *cv,
               double *n_sd, double *ests, double *work)
{
  int i, j, p, half, r, q, threep;
  double *m, *sd, *lam, *cut;
  int seed;

  /* calculate ns, len_plate */
  ns[0] = 0;
  for(p=0, r=0; p< *n_plates; p++, r+= *n_groups) {
    len_plate[p] = 0;
    for(i=0; i< *n_groups; i++) {
      len_plate[p] += n[i+r];
      ns[i+1+r] = ns[i+r] + n[i+r];
    }
  }

  /* split up workspace */
  m = work;
  sd = work + *n_plates * *n_groups;
  lam = sd + *n_plates;
  cut = lam + *n_plates * *n_plates;

  /* calculate means for each group in each plate */
  /* calculate SDs for cells alone group on each plate */
  /* calculate cut-off */
  for(p=0, r=0; p < *n_plates; p++, r+= *n_groups) {
    sd[p] = 0.0;
    for(i=0, q=0; i < *n_groups; i++, q+= *n_plates) {
      m[q + p] = 0.0;
      for(j=ns[i+r]; j<ns[i+1+r]; j++) {
        m[q + p] += y[j];
        if(!i) sd[p] += (y[j] * y[j]);
      }
      m[q + p] /= (double)(n[i+r]);
    }
    sd[p] = sqrt((sd[p] - (double)n[r] * m[p]*m[p])/(double)(n[r]-1));
    cut[p] = m[p] + *n_sd * sd[p];
  }

  /* find number of wells below cut-off for each group on each plate */
  /* then estimate lam */
  for(p=0, r=0; p < *n_plates; p++, r+= *n_groups) {
    for(i=0, q=0; i < *n_groups; i++, q+= *n_plates) {
      lam[q + p] = 0.0;
      for(j=ns[i+r]; j<ns[i+1+r]; j++) {
        if(y[j] < cut[p]) lam[q + p] += 1.0;
      }
      if(fabs(lam[q + p])< 0.01) lam[q + p] += 1.0;
      else if(fabs(lam[q+p] - (double)n[i+r]) < 0.01) lam[q + p] -= 1.0;
      lam[q + p] = - log(lam[q + p] / (double)n[i+r]) * c[p* *len_plate];
    }
  }
  /* sort lambdas within a group, then calculate median */
  for(i=0, q=0; i< *n_groups; i++, q+= *n_plates) {
    piksrt(*n_plates, lam + q);
    half = (*n_plates/2);
    if(fabs((double)(half*2 - *n_plates)) < 1e-6)
      ests[i] = (lam[q + half-1] + lam[q + half])/2.0;
    else ests[i] = lam[q + half];
  }

  /* find a_p:  median count among cells alone wells */
  /* place [mean(y in i,p) - a_p]/lam_i in lams */
  for(p=0, r=0, threep=*n_groups; p< *n_plates;
      p++, r+= *n_groups, threep += 3) {
    piksrt(n[r], y+ns[r]);
    half = (n[r]/2);
    if(fabs((double)(half*2 - n[r])) < 1e-6)
      ests[threep] = (y[ns[r] + half-1] + y[ns[r] + half])/2.0;
    else ests[threep] = y[ns[r] + half];

    for(i=0, q=0; i< *n_groups; i++, q+= *n_plates)
      lam[i +  r] = (m[q + p] - ests[threep]) / ests[i];
    piksrt(*n_groups, lam + r);
    half=(*n_groups/2);
    if(fabs((double)(half*2 - *n_groups)) < 1e-6)
      ests[threep+1] = (lam[r + half-1] +
                             lam[r + half])/2.0;
    else
      ests[threep+1] = lam[r + half];
    if(ests[threep+1] < 0.0) ests[threep+1] = 0.5;
    ests[threep+2] = ests[threep+1]/3.0;
  }

  if(fabs(*cv) > 1e-6) { /* randomize the starting point */
    seed = get_seed();
    for(i=0; i < *n_groups + 3 * *n_plates; i++) {
      ests[i] += gasdev(&seed)*fabs(ests[i])*(*cv);
      if(ests[i] < 0) ests[i] *= -1.0; /* none of the estimates should be negative */
    }
  }

}




/**********************************************************************
 *
 * piksrt
 *
 * This function was take from Numerical Recipes in C (Sec 8.1:
 * straight insertion and Shell's method).  It should suffice for my
 * very simple sorting needs.
 *
 * Input:
 *
 *     n = length of vector to sort
 *
 *   arr = pointer to array to be sorted into ascending order;
 *         on output, it contains the sorted array
 *
 **********************************************************************/

void piksrt(int n, double *arr)
{
  int i,j;
  double a;

  for(j=1; j<n; j++) {        /* pick out each element in turn */
    a=arr[j];
    i=j-1;
    while(i>=0 && arr[i]>a) {   /* look for the place to insert it */
      arr[i+1]=arr[i];
      i--;
    }
    arr[i+1]=a;                /* insert it */
  }
}



/****************************************
 *
 * get_seed
 *
 * returns a negative int, obtained from
 * the system clock
 *
 ****************************************/

int get_seed(void)
{
  time_t d1, d2;
  struct tm *e;

  d2 = time(&d1);
  e = localtime(&d1);
  return(-(((e->tm_mday*24 + e->tm_sec)*60 + e->tm_min)*60 + e->tm_hour));
}


/**********************************************************************
 *
 * gasdev
 *
 *   generates a pseudorandom normal(0,1) deviate
 *
 *   "idum" should be a int negative seed on first use
 *
 *   taken from numerical recipes in c
 *
 **********************************************************************/

double gasdev(int *idum)
{
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;
  double ran1();

  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}



/**********************************************************************
 *
 * ran1
 *
 *   generates a pseudorandom number uniformly from (0,1)
 *
 *   "idum" should be a int *negative* seed on first use
 *
 *   taken from numerical recipes in c
 *
 **********************************************************************/

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(int *idum)
{
  static int ix1,ix2,ix3;
  static double r[98];
  double temp;
  static int iff=0;
  int j;

  if (*idum < 0 || iff == 0) {
    iff=1;
    ix1=(IC1-(*idum)) % M1;
    ix1=(IA1*ix1+IC1) % M1;
    ix2=ix1 % M2;
    ix1=(IA1*ix1+IC1) % M1;
    ix3=ix1 % M3;
    for (j=1;j<=97;j++) {
      ix1=(IA1*ix1+IC1) % M1;
      ix2=(IA2*ix2+IC2) % M2;
      r[j]=(ix1+ix2*RM2)*RM1;
    }
    *idum=1;
  }
  ix1=(IA1*ix1+IC1) % M1;
  ix2=(IA2*ix2+IC2) % M2;
  ix3=(IA3*ix3+IC3) % M3;
  j=1 + ((97*ix3)/M3);
  if (j > 97 || j < 1) Rprintf("RAN1: This cannot happen.");
  temp=r[j];
  r[j]=(ix1+ix2*RM2)*RM1;
  return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3
