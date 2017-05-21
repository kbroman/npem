/**********************************************************************
 *
 * npem
 * version 0.50
 * 22 June 2003
 *
 * Karl W Broman
 *
 * npem_em.c : npem_em, npem_em_m, npem_em_e, npem_ll
 *
 * July, 1995 (revised 10/11/95, 2/3/96, 3/12/96, 4/16/96, 4/24/96,
 *                     8/18/2000, 6/22/2003)
 *
 * See the file "npem.h" for a thorough description.
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "npem.h"
#include <R.h>
#include <R_ext/PrtUtil.h>

/*********************************************************************
 *
 * npem_em
 *
 * This is the main function to perform the EM algorithm on the
 * normal/Poisson model described in npem.h
 *
 * Input: (all are pointers to either double or int)
 *
 *      y           the vector of scintillation counts y{pij}.
 *                  These should be in "lexicographic order,"
 *                  with the wells moving most quickly, then antigen
 *                  groups, and finally plates.  All plates should
 *                  contain at least one well for each group.
 *
 *      n_plates    the number of plates in the NPEM
 *
 *      c           the vector of cell densities c{pij}.
 *                  c{pij} gives the cell density corresponding to
 *                  y{pij}.  (Usually this is either a vector of 1's
 *                  or cells/well / 10^6 cells.)  [Note:  these can't
 *                  be 0.]
 *
 *      n_groups    the number of groups on each plate
 *
 *      n           n{pi} = number of wells for group i on plate p
 *                  (a vector of length n_plates x n_groups)
 *
 *        ns, len_plate, and len_group can be empty on input,
 *      they'll be calculated from n
 *
 *      ns          ns{pi} = partial sums n = (0,n[0],n[0]+n[1],...)
 *                  (one more than length of n)
 *
 *      len_plate   len_plate{p} = number of wells used on plate p
 *                  (vector of length n_plates) = sum_i n{pi}
 *
 *      len_group   len_group{i} = total number of wells for group i
 *                  (across all plates) = sum_p n{pi} (length n_groups)
 *
 *      ests        initial estimates (lambda's, (a, b, sigma)'s)
 *
 *      k           empty vector in which the calculated
 *                  E(k{pij}|y{pij}) will be placed (same length as y)
 *
 *      ksq         empty vector in which the calculated
 *                  E[k{pij}^2|y{pij}] will be placed (same length as y)
 *
 *      loglik      empty double in which the log likelihood at the MLE
 *                  will be placed
 *
 *      maxit       maximum number of iterations to perform
 *
 *      tol         tolerance value to determine when to stop
 *
 *      maxk        maximum value at which to perform sum when
 *                  calculating E(k{pij}|y{pij}) and E(k{pij}^2|y{pij})
 *
 *      prnt        0  don't print anything out
 *                  1  print the increase in the log lik. at each step
 *                  2  print the current estimates as well as the
 *                     increase in the log likelihood at each step
 *
 *      intwork    workspace of length n_groups
 *
 *      doublework  workspace of length 3(n_groups) + 8*n_plates
 *
 *      use_order_constraint   0  no constraints
 *                             1  if you want lambda_1 <= lambda_i
 *                             2  especially for the NIH 1394 data
 *                             3  lambda_i <= lambda_(i + n_groups/2)
 *                                (especially for the D. Koelle data)
 *
 * Output: void
 *
 **********************************************************************/

void npem_em(double *y, int *n_plates, double *c, int *n_groups,
            int *n, int *ns, int *len_plate, int *len_group,
            double *ests, double *k, double *ksq, double *loglik,
            int *maxit, double *tol, int *maxk, int *prnt,
            int *intwork, double *doublework,
            int *use_order_constraint)
{
  int i, j, p, flag, n_param;
  double oldll, *oldests;
  int *ns_cur, *n_cur;

  /* calculate ns, len_plate, len_group */
  ns[0] = 0;
  for(i=0; i < *n_groups; i++) len_group[i] = 0;
  ns_cur = ns; n_cur = n;
  for(p=0; p < *n_plates; p++) {
    len_plate[p] = 0;
    for(i=0; i< *n_groups; i++, ns_cur++, n_cur++) {
      len_plate[p] += *n_cur;
      len_group[i] += *n_cur;
      *(ns_cur+1) = *ns_cur + *n_cur;
    }
  }

  oldests = doublework;
  n_param = *n_plates * 3 + *n_groups;

  /* make sure the lambdas are not strictly 0 */
  for(i=0; i< n_param; i++) if(ests[i] < 1e-50) ests[i]=1e-50;

  /* calculate initial loglik */
  if(*prnt)
    npem_ll(y, n_plates, c, n_groups, n, ns, len_plate, len_group,
           ests, loglik, maxk);

  /* E step */
  npem_em_e(y, n_plates, c, n_groups, n, ns, ests, k, ksq, maxk);

  for(i=0; i<*maxit; i++) {
    for(j=0; j< n_param; j++) oldests[j] = ests[j];

    /* one EM step */
    npem_em_m(y, n_plates, c, n_groups, n, ns, len_plate, ests,
             k, ksq, intwork, doublework, use_order_constraint);
    npem_em_e(y, n_plates, c, n_groups, n, ns, ests, k, ksq, maxk);

    /* check for convergence */
    for(j=0, flag=1; j<n_param; j++) {
      if(fabs(ests[j]-oldests[j]) > *tol) {
        flag=0; break;
      }
    }
    if(flag) break;

    /* print out stuff if desired */
    if(*prnt == 2) for(j=0; j<n_param; j++)
      Rprintf("  %6.2lf", ests[j]);
    if(*prnt) {
      oldll = *loglik;
      npem_ll(y, n_plates, c, n_groups, n, ns, len_plate, len_group,
             ests, loglik, maxk);
      Rprintf("  %8.6lf\n", *loglik - oldll);
    }
  }
  *maxit = i;
  /* calculate final likelihood */
  npem_ll(y, n_plates, c, n_groups, n, ns, len_plate, len_group,
         ests, loglik, maxk);

  Rprintf("\tnumber of iterations: %ld\n", i+1);
}

/**********************************************************************
 * npem_em_e
 *
 *  This function performs the E step of the EM algorithm for
 *  the above-described model.
 *
 * Input: (all are pointers to either double or int)
 *
 *      y, n_plates, c, n_groups, n, ns, ests, k, ksq, maxk
 *
 *      (all just like in npem_em)
 *
 * Output: void
 *
 **********************************************************************/

void npem_em_e(double *y, int *n_plates, double *c, int *n_groups,
              int *n, int *ns, double *ests, double *k, double *ksq,
              int *maxk)
{
  int i, j, p, r;
  double f, g, h, du;
  int *ns_cur;
  int absig_loc;

  /* calculate E(k_ij|y_ij) and E[(k_ij)^2|y_ij] */

  ns_cur = ns;
  absig_loc = *n_groups;

  /* p goes over the plates */
  for(p=0; p < *n_plates; p++, absig_loc += 3) {

    /* i goes over lambda groups within a plate */
    for(i=0; i< *n_groups; i++, ns_cur++) {

      /* j goes over the y within a lambda group */
      for(j= *ns_cur; j< *(ns_cur + 1); j++) {
        k[j] = 0.0;
        ksq[j] = 0.0;

        h = (y[j]- ests[absig_loc]) / ests[absig_loc+2];
        g = exp( -0.5 * h * h );

        for(f=0.0, r=1; r < (*maxk); r++) {

          h = (y[j] - ests[absig_loc] - ests[absig_loc+1] *
               ((double)r)) / ests[absig_loc+2];

          f = ( -0.5 * h * h + ((double)r) * log(ests[i]*c[j]) );
          du = gammaln((double)r);
          k[j] += exp(f - du);
          ksq[j] += exp(f + log((double)r) - du);
      du = gammaln((double)(r+1));
          g += exp(f - du);
        }
        if(g < 1e-50) {
          k[j]=1e-50;              /* this hopefully fixes a problem */
          ksq[j]=1e-50;            /*   which occurs when g=0        */
        }
        else {
          k[j] /= g;
          ksq[j] /= g;
        }
      }
    }
  }
}


/**********************************************************************
 *
 * npem_em_m
 *
 * This function performs the M-step of the EM algorithm for the above-
 * described model.
 *
 * Input: (all are pointers to either double or int)
 *
 *      y, n_plates, c, n_groups, n, ns, len_plate, ests,
 *      k, ksq, intwork, doublework, use_order_constraint
 *
 *      (all just like in npem_em)
 *
 * Output: void
 *
 **********************************************************************/

void npem_em_m(double *y, int *n_plates, double *c, int *n_groups,
              int *n, int *ns, int *len_plate, double *ests,
              double *k, double *ksq, int *intwork,
              double *doublework, int *use_order_constraint)
{
  int i, j, i2, p, *which, n_param;
  double *sy, *ssk, *ssy, *sky, *sk;
  double check, *sk_groups, *sc_groups, sumc;
  double lpd;
  int *ns_cur;
  int absig_loc;

  /* setup work space */
  n_param = *n_groups + 3* *n_plates;
  sk_groups = doublework + *n_plates*3 + *n_groups;
  sc_groups = sk_groups + *n_groups;
  sk = sc_groups + *n_groups;
  sy = sk + *n_plates;
  sky = sy + *n_plates;
  ssk = sky + *n_plates;
  ssy = ssk + *n_plates;
  which = intwork;

  /* obtain sufficient statistics */
  for(i=0; i<*n_groups; i++) {
    sk_groups[i] = 0.0;
    sc_groups[i] = 0.0;
  }

  ns_cur=ns;
  for(p=0; p<*n_plates; p++) { /* p varies across plates */
    sy[p]=0.0;
    sky[p]=0.0;
    ssk[p]=0.0;
    ssy[p]=0.0;
    sk[p]=0.0;
    /* i varies across groups */
    for(i=0; i<*n_groups; i++,ns_cur++) {
      for(j= *ns_cur; j< *(ns_cur+1); j++) {
        sy[p] += y[j];
        sky[p] += (y[j]*k[j]);
        ssk[p] += ksq[j];
        ssy[p] += (y[j]*y[j]);
        sk[p] += k[j];
        sk_groups[i] += k[j]; /* sum of k's within a group */
        sc_groups[i] += c[j]; /* sum of concentrations within a group */
      }
    }
  }

  /* lambda's obtained as means of k's (weighted by concentrations) */
  for(i=0; i< *n_groups; i++)
    ests[i] = sk_groups[i] / sc_groups[i];

  if(*use_order_constraint==3) {
    /* check and fix order constraints */
    /* this is especially for the David Koelle data */

    i = *n_groups/2;
    for(i2=0; i2<i; i2++) {
      if(ests[i2+i]<ests[i2]) {
        ests[i2]=sk_groups[i2] + sk_groups[i+i2];
        sumc = sc_groups[i2] + sc_groups[i+i2];
        ests[i2] /= sumc;
        ests[i+i2] = ests[i2];
      }
    }
  }
  else if(*use_order_constraint==2) {
    /* check and fix order constraints */
    /* this is especially for the NIH 1394 data */

    for(i2=0; i2<16; i2 += 4) {
      for(i=0; i<3;i++) {
        which[i]=0; check=9e+99;
        for(j=1; j<4; j++) {
          if(ests[j+i2]<ests[i2] && ests[j+i2]<check) {
            check=ests[j+i2];
            which[i]=j;
          }
        }
        if(!which[i]) break;
        ests[i2]=sk_groups[i2]; sumc = sc_groups[i2];
        for(j=0; j<i+1; j++) {
          ests[i2] += sk_groups[which[j]+i2];
          sumc += sc_groups[which[j]+i2];
        }
        ests[i2] /= sumc;
        for(j=0; j<i+1; j++) ests[which[j]+i2] = ests[i2];
      }
    }
  }
  else if(*use_order_constraint) {
    /* check and fix order constraints */
    for(i=0; i<*n_groups-1;i++) {
      which[i]=0; check=9e+99;
      for(j=1; j<*n_groups; j++) {
        if(ests[j]<ests[0] && ests[j]<check) {
          check=ests[j];
          which[i]=j;
        }
      }
      if(!which[i]) break;
      ests[0]=sk_groups[0]; sumc = sc_groups[0];
      for(j=0; j<i+1; j++) {
        ests[0] += sk_groups[which[j]];
        sumc += sc_groups[which[j]];
      }
      ests[0] /= sumc;
      for(j=0; j<i+1; j++) ests[which[j]] = ests[0];
    }
  }

  /* lower threshold for the lambdas, to ensure that loglik can be calc. */
  for(i=0; i< *n_groups; i++) if(ests[i] < 1e-50) ests[i] = 1e-50;

  /* LS using E(k|y) and E(k^2|y) */
  absig_loc = *n_groups;
  for(p=0; p<*n_plates; p++, absig_loc += 3) {

    lpd = (double) len_plate[p];
    ests[absig_loc+1] =
      (sky[p] - sy[p]*sk[p]/lpd)/
        (ssk[p] - sk[p]*sk[p]/lpd);
    ests[absig_loc] =
      (sy[p] - ests[absig_loc+1]*sk[p])/lpd;
    ests[absig_loc+2] =
      sqrt((ssy[p]-sy[p]*sy[p]/lpd-ests[absig_loc+1]*
            (sky[p]-sy[p]*sk[p]/lpd))/lpd);
  }
  /* lower threshold for the a, b, sigma */
  for(i= *n_groups; i< n_param; i++) if(ests[i] < 1e-50) ests[i] = 1e-50;

}

/***********************************************************************
 *
 * npem_ll
 *
 * This function calculates the log likelihood for a set of parameter
 * values given a set of y data for the above-described normal/
 * Poisson model.
 *
 * Input: (all are pointers to either double or int)
 *
 *      y, n_plates, c, n_groups, n, ns, len_plate, len_group,
 *      ests, loglik, maxk
 *
 *      (all just like in npem_em)
 *
 * Output: void
 *
 ***********************************************************************/

void npem_ll(double *y, int *n_plates, double *c, int *n_groups,
            int *n, int *ns, int *len_plate, int *len_group,
            double *ests, double *loglik, int *maxk)
{
  int i, j, r, p;
  double f, g, du;
  int *ns_cur;
  int absig_loc;

  /* calculate log likelihood */
  *loglik = 0.0;

  ns_cur = ns;
  absig_loc = *n_groups;
  /* p goes through the plates */
  for(p=0; p<*n_plates; p++, absig_loc += 3) {

    /* i goes through the lambda groups */
    for(i=0; i<*n_groups; i++, ns_cur++) {

      /* j goes through the y within a lambda group */
      for(j= *ns_cur; j< *(ns_cur+1); j++) {
        f = 0.0;

        for(r=0; r < (*maxk); r++) {
          g = (y[j] - ests[absig_loc] - ests[absig_loc+1]
               * ((double)r)) / ests[absig_loc+2];
      du = gammaln((double)(r+1));
          f += exp(-0.5*g*g + ((double)r)*log(ests[i]*c[j]) - du);
        }
        if(f <= 1e-50) f = 1e-50;
        (*loglik) += (log(f) - ests[i]*c[j]);
      }
      (*loglik) -= ( (double)n[i] * log(ests[absig_loc+2]) );
    }
  }
}



/* returns the log of the gamma function */
/* (taken from Numerical Recipes in C) */
double gammaln(double xx)
{
  double x, y, tmp, ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2, -0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for(j=0; j<=5; j++) ser += cof[j]/++y;
  return(-tmp+log(2.5066282746310005*ser/x));
}
