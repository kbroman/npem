/**********************************************************************
 *
 * npem
 * version 0.50
 * 22 June 2003
 *
 * Karl W Broman
 *
 * npem_sem.c : npem_sem, npem_full_inf, npem_get_rates, npem_get_var
 *
 * July, 1995 (revised 10/11/95, 2/3/96, 3/12/96, 3/13/96,
 *                      4/16/96, 4/24/96, 4/26/96, 4/30/96,
 *                      5/3/96, 8/18/2000, 6/22/2003)
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
 * npem_sem
 *
 * This is the main function to perform the SEM algorithm on the
 * normal/Poisson model described in npem.h, to get an estimated
 * var-cov matrix.  Only the full-data information matrix and the
 * rate matrix are calculated.
 *
 * Input: (all are pointers to either double or int)
 *
 *      y           the vector of scintillation counts y{pij}.
 *                  These should be in "lexicographic order,"
 *                  with the wells moving most quickly, then groups,
 *                  and finally plates.
 *
 *      n_plates    the number of plates in the NPEM
 *
 *      c           the vector of cell densities c{pij}.
 *                  c{pij} gives the cell density corresponding to
 *                  y{pij}.  (They should give cells/well / 10^6 cells.)
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
 *      n_se        the number of se's to obtain
 *                  (Generally this is n_groups + 3*n_plates, but it
 *                  could be just n_groups, if we wish to ignore the
 *                  (a,b,sigma), to make things easier.)
 *                  Actually n_se can be only either n_groups or
 *                  n_groups + 3*n_plates, and nothing else.
 *
 *      mle         the mle, as obtained from npem_em, in the form
 *                  (lambda's, (a, b, sigma)'s)
 *
 *      start       a vector similar to "mle," but some short distance
 *                  away, at which the SEM algorithm begins.  (We
 *                  monitor the rate at which the EM algorithm brings
 *                  us to the mle, in order to get the rate matrix.
 *
 *      k           empty vector in which the calculated
 *                  E(k{pij}|y{pij}) will be placed (same length as y)
 *
 *      ksq         empty vector in which the calculated
 *                  E[k{pij}^2|y{pij}] will be placed (same length as y)
 *
 *      infor       empty square matrix of doubles, with dimension
 *                  (n_se x n_se) in which the full-data information
 *                  matrix will be placed
 *
 *      rates       empty square matrix of doubles, the same size as
 *                  "infor," in which the rate matrix ("DM" in Meng
 *                  and Rubin's notation) will be placed
 *
 *      var         empty square matrix of doubles, the same size as
 *                  "infor," in which the estimated var-cov matrix will
 *                  be placed.
 *
 *      se          empty vector of length n_se, in which the estimated
 *                  standard errors will be placed.
 *
 *      get_var     0 don't calculate var-cov matrix and se's
 *                  1 do calculate var-cov matrix and se's
 *
 *      tol         tolerance value to determine when to stop
 *
 *      maxk        maximum value at which to perform sum when
 *                  calculating E(k{pij}|y{pij}) and E(k{pij}^2|y{pij})
 *
 *      prnt        0  don't print anything out
 *                  1  print the rate matrix and which elements have
 *                     converged at each step
 *
 *      intwork    workspace of length 10*(n_param*n_param + n_param)
 *
 *      doublework  workspace of length 20*(n_param*n_param + n_param)
 *
 *      use_order_constraint   1 if you want lambda_1 <= lambda_i
 *
 *      maxit       int giving maximum number of iterations before
 *                  exit in the npem_get_rates function
 *
 * Output: void
 *
 **********************************************************************/

void npem_sem(double *y, int *n_plates, double *c, int *n_groups,
             int *n, int *ns, int *len_plate, int *len_group,
             int *n_se, double *mle, double *start,
             double *k, double *ksq, double *infor, double *rates,
             double *var, double *se, int *get_var,
             double *tol, int *maxk, int *prnt, int *intwork,
             double *doublework, int *use_order_constraint,
             int *maxit)
{
  int i, p;
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

  npem_full_inf(y,n_plates,c,n_groups,n,ns,len_plate, len_group,
               n_se,mle,k,ksq,infor,doublework);
  npem_get_rates(y,n_plates,c,n_groups,n,ns,len_plate, len_group,
                n_se,mle,start,k,ksq,rates,
                tol,maxk,prnt,intwork,doublework,
                use_order_constraint, maxit);
  if(*get_var)
    npem_get_var(infor, rates, n_se, var, se, intwork, doublework);
}

/**********************************************************************
 *
 * npem_full_inf
 *
 * The following calculates the full-data information matrix.
 *
 * Input: (all are pointers to either double or int)
 *
 *      y, n_plates, c, n_groups, n, ns, len_plate, len_group,
 *      n_se, mle, k, ksq, infor, doublework
 *
 *      (all just like in npem_sem)
 *
 * Output: void
 *
 **********************************************************************/

void npem_full_inf(double *y, int *n_plates, double *c, int *n_groups,
                  int *n, int *ns, int *len_plate, int *len_group,
                  int *n_se, double *mle, double *k,
                  double *ksq, double *infor, double *doublework)
{

  int i, j, p;              /* indices */
  int n_param;              /* # of parameters */
  double ym=0.0;             /* mean of y's within a plate */
  double s1, s2, s5;         /* various sums within a plate */
  double sigmasq;            /* sigma squared */
  double *sc_groups;         /* sum(c) within a lambda group */
  int end_of_prev_plate, end_of_cur_plate;  /* used in loop over plate */
  int *ns_cur;
  int absig_loc;

  /* set up workspace, etc. */
  n_param = *n_groups + 3 * *n_plates;
  sc_groups=doublework;

  /* calculate sum(c) within a "lambda group" */
  for(i=0; i< *n_groups; i++) sc_groups[i]=0.0;
  ns_cur = ns;
  for(p=0; p < *n_plates; p++)
    for(i=0; i< *n_groups; i++, ns_cur++)
      for(j=*ns_cur; j < *(ns_cur+1); j++)
        sc_groups[i] += c[j];

  /* fill information matrix with 0's */
  for(i=0; i< *n_se * *n_se; i++) infor[i] = 0.0;

  /* ensure that mle's are > 1e-50 */
  for(i=0; i< n_param; i++)
    if(mle[i] < 1e-50) mle[i] = 1e-50;

  /* calculate first n_groups diag elements */
  for(i=0; i< *n_groups; i++)
    infor[i+ *n_se*i] = sc_groups[i]/mle[i];

  if(*n_se == n_param) {

    /* calculate the rest of the information matrix */
    end_of_prev_plate=0; end_of_cur_plate=0;
    absig_loc = *n_groups;
    for(p=0; p < *n_plates; p++, absig_loc += 3) {
      end_of_cur_plate += len_plate[p];

      sigmasq = mle[absig_loc+2]*mle[absig_loc+2];

      infor[absig_loc + n_param * absig_loc] =
        (double)len_plate[p] / sigmasq;

      for(i=end_of_prev_plate; i<end_of_cur_plate; i++) ym += y[i];
      ym /= (double)len_plate[p];

      /* calculate various sums with a plate */
      s1 = 0.0; s2 = 0.0; s5 = 0.0;
      for(i=end_of_prev_plate;i<end_of_cur_plate; i++) {
        s1 += k[i];
        s2 += ksq[i];
        s5 += ((y[i] - ym)*(y[i]-ym-mle[absig_loc+1]*k[i]));
      }

      infor[(absig_loc+1)+n_param*(absig_loc+1)] = s2 / sigmasq;

      infor[(absig_loc+2)+n_param*(absig_loc+2)] =
        (3.0*s5/sigmasq - ((double)len_plate[p]))/sigmasq;

      infor[(absig_loc)+(n_param)*(absig_loc+1)] =
        infor[(absig_loc+1)+(n_param)*(absig_loc)] = s1 / sigmasq;

      end_of_prev_plate += len_plate[p];
    }
  }
}

/**********************************************************************
 *
 * npem_get_rates
 *
 * The following calculates a rate matrix as part of the SEM algorithm
 *
 * Input: (all are pointers to either double or int)
 *
 *      y, n_plates, c, n_groups, n, ns, len_plate, len_group,
 *      n_se, mle, start, k, ksq,
 *      rates, tol, maxk, prnt, intwork, doublework,
 *      use_order_constraint, maxit
 *
 *      (all just like in npem_sem)
 *
 * Output: void
 *
 **********************************************************************/

void npem_get_rates(double *y, int *n_plates, double *c, int *n_groups,
                   int *n, int *ns, int *len_plate, int *len_group,
                   int *n_se, double *mle, double *start, double *k,
                   double *ksq, double *rates, double *tol, int *maxk,
                   int *prnt, int *intwork, double *doublework,
                   int *use_order_constraint, int *maxit)
{
  int i, j, *num_converge, *flags, flag2, *intwork2;
  double *next, *dummy, *oldrates, *doublework2;
  int n_param, st;
  int n_se_sq;
  int cur_loc;

  /* number of parameters */
  n_param = *n_groups + 3 * *n_plates;

  /* divvy up workspace */
  num_converge = intwork;
  flags = num_converge + n_param;
  intwork2 = flags + n_param * n_param;
  next = doublework;
  dummy = next + n_param;
  oldrates = dummy + n_param;
  doublework2 = oldrates + n_param*n_param;

  for(i=0; i< n_param; i++) {
    next[i] = start[i];
    num_converge[i]=0;
  }
  for(i=0; i< *n_se * *n_se; i++) {
    rates[i] = -50000.0;
    flags[i]=0;
  }

  n_se_sq = *n_se * *n_se;

  for(st=0; st< *maxit; st++) {
    for(i=0; i< n_se_sq; i++) oldrates[i] = rates[i];

    /* take one EM step */
    npem_em_e(y, n_plates, c, n_groups, n, ns, next, k, ksq, maxk);
    npem_em_m(y, n_plates, c, n_groups, n, ns, len_plate,
             next, k, ksq, intwork2, doublework2, use_order_constraint);

    flag2=0;    /* this flag indicates that not all of the */
                /* rates have converged */

    for(i=0; i< *n_se; i++) {

      /* calculate row i of the matrix, if not all of its elements converged */
      if(num_converge[i] < *n_se) {

        /* form mle vector with ith element replaced by the current step */
        for(j=0; j< n_param; j++) dummy[j]=mle[j];
        dummy[i] = next[i];

        /* take one EM step */
        npem_em_e(y, n_plates, c, n_groups, n, ns, dummy, k, ksq, maxk);
        npem_em_m(y, n_plates, c, n_groups, n, ns, len_plate,
                 dummy, k, ksq, intwork2, doublework2, use_order_constraint);

        for(j=0; j< *n_se; j++) {
          cur_loc = i + *n_se * j;

          /* calculate rate, if this one hasn't converged */
          if(!flags[cur_loc]) {
            if(fabs(next[i]-mle[i])<1e-50)
              rates[cur_loc] = 1e+20;
            else
              rates[cur_loc] = (dummy[j]-mle[j])/(next[i]-mle[i]);

            /* check for convergence */
            if(fabs(rates[cur_loc] -
                    oldrates[cur_loc])<(*tol)) {
              flags[cur_loc] = 1;
              num_converge[i]++;
            }
          }
        }
        flag2=1;
      }
    }
    if(!flag2) return;              /* all rates converged */

    if(*prnt) {
      for(i=0; i< *n_se; i++) {
        for(j=0; j< *n_se; j++)
          Rprintf("%5.2lf %1ld ",rates[i+n_param*j],
                  flags[i+*n_se*j]);
        Rprintf("\n");
      }
      Rprintf("\n");
    }
  }
  *maxit = st;
}

/**********************************************************************
 *
 * npem_get_var
 *
 *      This function takes the full-data information matrix and the
 *   matrix of rates, calculated above, and obtains the estimated
 *   variance-covariance matrix.
 *
 *   INPUT:
 *
 *     infor = matrix of size (n_param x n_param), obtained from
 *             "npem_full_infor"
 *
 *     rates = matrix of size (n_param x n_param), obtained from
 *             "npem_get_rates"
 *
 *     n_param = number of parameters (n_se in "npem_sem")
 *
 *     var = matrix of size (n_param x n_param), in which the estimated
 *           variance-covariance matrix will be placed
 *
 *     se = vector of length n_param, in which the estimated SEs will
 *          be placed
 *
 *     intwork = vector of length 3 * n_param to use as workspace
 *
 *     doublework = vector of length (4 * n_param*n_param + n_param)
 *                  to use as workspace
 *
 **********************************************************************/

void npem_get_var(double *infor, double *rates, int *n_param,
                 double *var, double *se, int *intwork,
                 double *doublework)
{
  int i, j, k, flag1, psq, errorflag;
  int *bad, *good, n_bad, n_good, *inttemp;
  double *temp1, *temp2, *temp3, *temp4, *temp5, var_diff=0.0;

  psq = *n_param * *n_param;

  /* divvy up workspace */
  bad = intwork;
  good = bad + *n_param;
  inttemp = good + *n_param;
  temp1 = doublework;
  temp2 = temp1 + psq;
  temp3 = temp2 + psq;
  temp4 = temp3 + psq;
  temp5 = temp4 + psq;

  /* place infor in temp2 */
  for(i=0; i<psq; i++) temp2[i] = infor[i];

  /* invert temp2, place result in temp1 */
  invert_matrix(temp2, temp1, n_param, inttemp, temp5, &errorflag);
  if(errorflag) {
    Rprintf("Error when inverting full-data information matrix\n");
    return;
  }

  /* find rows which appear to be "bad", containing really high rates */
  n_bad = 0; n_good = 0;
  for(i=0; i< *n_param; i++) {
    for(j=0, flag1=0; j< *n_param; j++, flag1=0) {
      if(fabs(rates[i + *n_param * j]) >= 1e+5) {
        flag1 = 1;
        bad[n_bad] = i;
        n_bad ++;
        break;
      }
    }
    if(!flag1) {
      good[n_good] = i;
      n_good++;
    }
  }

  if(!n_bad) {  /* if no bad rows */

    /* place I - rates in temp3 */
    for(i=0; i< psq; i++) temp3[i] = -rates[i];
    for(i=0; i< *n_param; i++) temp3[i + *n_param * i] += 1.0;

    /* invert temp3 and place in temp2 */
    invert_matrix(temp3, temp2, n_param, inttemp, temp5, &errorflag);
    if(errorflag) {
      Rprintf("Error when inverting (I-rates) matrix\n");
      return;
    }

    /* let temp4 = rates * temp2 */
    for(i=0; i< *n_param; i++)
      for(j=0; j< *n_param; j++)
        for(k=0, temp4[i+ *n_param * j] = 0.0; k < *n_param; k++)
          temp4[i+ *n_param*j] += rates[i+ *n_param*k] * temp2[k+ *n_param*j];

    /* let temp4 = I + temp4 */
    for(i=0; i< *n_param; i++) temp4[i+ *n_param*i] += 1.0;

    /* let var = temp1 * temp4 */
    for(i=0; i<*n_param; i++)
      for(j=0; j<*n_param; j++)
        for(k=0, var[i+ *n_param * j] = 0.0; k < *n_param; k++)
          var[i+ *n_param*j] += temp1[i+ *n_param*k] * temp4[k+ *n_param*j];

  }

  else {   /* some "bad" rows */

/**********************************************************************
       let X     = inverse of "infor" (now in temp1)
       let G1    = "bad" rows and columns of X
       let G2    = "bad" rows, "good" columns of X
       let G3    = "good" rows and columns of X
       let Rstar = "good" rows and columns of "rates"

       Then we wish to calculate:
               V = X + (G3 - G2' inverse(G1) G2) Rstar inverse(I-Rstar)
 **********************************************************************/

    /* let temp2 = I - Rstar */
    for(i=0; i<n_good; i++) {
      for(j=0; j<n_good; j++)
        temp2[i+ n_good * j] = -rates[good[i] + *n_param * good[j]];
      temp2[i+ n_good * i] += 1.0;
    }

    /* invert temp2 and place in temp3 */
    invert_matrix(temp2, temp3, &n_good, inttemp, temp5, &errorflag);
    if(errorflag) {
      Rprintf("Error when inverting (I - rates*) matrix\n");
      return;
    }

    /* let temp2 = Rstar * temp3 */
    for(i=0; i<n_good; i++)
      for(j=0; j<n_good; j++)
        for(k=0, temp2[i+n_good*j] = 0.0; k<n_good; k++)
          temp2[i+n_good*j] +=
            rates[good[i]+ *n_param*good[k]] * temp3[k+n_good*j];

    /* let temp3 = G1 (bad rows and columns of temp1) */
    for(i=0; i<n_bad; i++)
      for(j=0; j<n_bad; j++)
        temp3[i+n_bad*j] = temp1[bad[i] + *n_param * bad[j]];

    /* invert temp3 and place in temp4 */
    invert_matrix(temp3, temp4, &n_bad, inttemp, temp5, &errorflag);
    if(errorflag) {
      Rprintf("Error when inverting part of information matrix\n");
      return;
    }

    /* let temp3 = temp4 * G2 (bad rows, good columns of temp1) */
    for(i=0; i<n_bad; i++)
      for(j=0; j<n_good; j++)
        for(k=0, temp3[i+n_bad*j]=0.0; k<n_bad; k++)
          temp3[i+n_bad*j] +=
            temp4[i+n_bad*k] * temp1[bad[k]+ *n_param*good[j]];

    /* let temp4 = G2' * temp3 */
    for(i=0; i<n_good; i++)
      for(j=0; j<n_bad; j++)
        for(k=0, temp4[i+n_good*j]=0.0; k<n_good; k++)
          temp4[i+n_good*j] +=
            temp1[bad[k]+ *n_param*good[i]]* temp3[k+n_bad*j];

    /* let temp3 = (G3 - temp4) * temp2 */
    for(i=0; i<n_good; i++)
      for(j=0; j<n_good; j++)
        for(k=0, temp3[i+n_good*j]=0.0; k<n_good; k++)
          temp3[i+n_good*j] =
            (temp1[good[i]+ *n_param*good[k]] - temp4[i+n_good*k]) *
              temp2[k+n_good*j];

    /* let var = temp1 */
    for(i=0; i<psq; i++) var[i] = temp1[i];

    /* add temp3 to the good rows and columns of var */
    for(i=0; i<n_good; i++)
      for(j=0; j<n_good; j++)
        var[good[i]+ *n_param*good[j]] += temp3[i+n_good*j];

    /* Done! */
  }

  /* calculate SE's */
  for(i=0, flag1=0; i< *n_param; i++)
    if(var[i+ *n_param*i] >= 0.0) se[i] = sqrt(var[i+ *n_param*i]);
    else {
      flag1=1;
      se[i] = -1.0;
    }
  if(!flag1) {
    /* check symmetry of estimated variance-covariance matrix */
    /* then symmetrize the matrix */
    for(i=0, flag1=0; i< *n_param; i++) {
      for(j=0; j<i; j++) {
        if((fabs(var[i + *n_param *j] - var[j + *n_param*i]) > 0.01)
           && (fabs(var[i + *n_param *j] - var[j + *n_param*i]) >
               sqrt(var[i + *n_param * i] * var[j + *n_param *j]))) {
          flag1=1;
          if(fabs(var[i+*n_param*j] - var[j+*n_param*i]) > var_diff)
            var_diff = fabs(var[i+*n_param*j] - var[j+*n_param*i]);
        }

        var[i + *n_param *j] = var[j + *n_param * i] =
          (var[i + *n_param *j] + var[j + *n_param * i])/2.0;
      }
    }

    if(flag1) {
      Rprintf("SEM algorithm has not properly converged.\n");
      Rprintf("(The estimated var-cov matrix is not symmetric.  ");
      Rprintf("The largest deviation is %0.3f)\n\n", var_diff);
    }
  }
  else {
    Rprintf("Some estimated variances are negative.\n");

    /* check symmetry of estimated variance-covariance matrix */
    /* then symmetrize the matrix */
    for(i=0, flag1=0; i< *n_param; i++) {
      for(j=0; j<i; j++) {
        if(fabs(var[i + *n_param *j] - var[j + *n_param*i]) > 0.01) {
          flag1=1;
          if(fabs(var[i+*n_param*j] - var[j+*n_param*i]) > var_diff)
            var_diff = fabs(var[i+*n_param*j] - var[j+*n_param*i]);
        }

        var[i + *n_param *j] = var[j + *n_param * i] =
          (var[i + *n_param *j] + var[j + *n_param * i])/2.0;
      }
    }

    if(flag1) {
      Rprintf("SEM algorithm has not properly converged.\n");
      Rprintf("(The estimated var-cov matrix is not symmetric.  ");
      Rprintf("The largest deviation is %0.3f)\n\n", var_diff);
    }
  }
}
