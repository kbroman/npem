/**********************************************************************
 *
 * npem
 * version 0.50
 * 2021-07-14
 *
 * Karl W Broman
 *
 * npem_start.c : npem_start, piksrt
 *
 * 3/13/96 (revised 4/25/96, 8/20/2000, 6/22/2003, 2021-07-14)
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
    R_rsort(lam + q, *n_plates);
    half = (*n_plates/2);
    if(fabs((double)(half*2 - *n_plates)) < 1e-6)
      ests[i] = (lam[q + half-1] + lam[q + half])/2.0;
    else ests[i] = lam[q + half];
  }

  /* find a_p:  median count among cells alone wells */
  /* place [mean(y in i,p) - a_p]/lam_i in lams */
  for(p=0, r=0, threep=*n_groups; p< *n_plates;
      p++, r+= *n_groups, threep += 3) {
    R_rsort(y+ns[r], n[r]);
    half = (n[r]/2);
    if(fabs((double)(half*2 - n[r])) < 1e-6)
      ests[threep] = (y[ns[r] + half-1] + y[ns[r] + half])/2.0;
    else ests[threep] = y[ns[r] + half];

    for(i=0, q=0; i< *n_groups; i++, q+= *n_plates)
      lam[i +  r] = (m[q + p] - ests[threep]) / ests[i];
    R_rsort(lam + r, *n_groups);
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
    GetRNGstate();

    for(i=0; i < *n_groups + 3 * *n_plates; i++) {
      ests[i] += norm_rand()*fabs(ests[i])*(*cv);
      if(ests[i] < 0) ests[i] *= -1.0; /* none of the estimates should be negative */
    }

    PutRNGstate();
  }

}
