/*********************************************************************
 *
 * npem
 * version 0.50
 * 22 August 2000
 *
 * Karl W Broman
 *
 * npem.h  (header file for npem_em.c, npem_sem.c, npem_start.c and
 *                            npem_matrix.c)
 *
 * July, 1995 (revised 10/11/95, 2/3/96, 3/12/96, 3/13/96, 4/26/96,
 *                      4/30/96, 5/3/96, 8/20/2000)
 *
 * These functions are written to perform the EM and SEM algorithms
 * on a limiting dilution assay (or even a single assay).  The data
 * come in the form of transformed scintillation counts, y{pij},
 * where p=plate, i=group, j=well.
 *
 * Well j of group i of plate p is assumed to contain k{pij}
 * responding cells, and gives the transformed scint. count y{pij}.
 * We assume that k{pij} ~ Poisson(lambda{i} c{pij}), where c{pij}
 * is a known concentration of cells.  y{pij}|k{pij} is assumed to
 * be normal(a{p} + b{p} k{pij}, sigma{p}^2).  We seek maximum
 * likelihood estimates of the lambda{i} and the (a{p},b{p},sigma{p})
 * using the EM algorithm.  Estimated SEs are obtained via the SEM
 * algorithm.
 *
 *
 * References:
 *
 *   Dempster A, Laird N, and Rubin D (1977) Maximum likelihood
 *   estimation from incomplete data via the EM algorithm.  Journal
 *   of the Royal Statistical Society, Series B. 39: 1-38.
 *
 *   Meng X-L and Rubin DB (1991) Using EM to obtain asymptotic
 *   variance-covariance matrices: the SEM algorithm.  Journal of the
 *   American Statistical Association. 86: 899-909.
 *
 *   Broman K, Speed T, and Tigges M (1996) Estimation of antigen-
 *   responsive T cell frequencies in PBMC from human subjects.
 *   J Immunol Meth 198:119-132.
 *
 *   Broman K, Speed T, and Tigges M (1998) Estimation of antigen-
 *   responsive T cell frequencies in PBMC from human subjects.
 *   Stat Sci 13:4-8
 *
 *********************************************************************
 *
 * functions:
 *
 *   npem_em.c:    npem_em        main EM program
 *                 npem_em_e      E step of EM program
 *                 npem_em_m      M step of EM program
 *                 npem_ll        calculates log likelihood
 *                 gammaln        calculates the log of the gamma
 *                                function (used to calculate n!)
 *
 *   npem_sem.c:   npem_sem       main SEM program
 *                 npem_full_inf  calculates full-data information
 *                                matrix
 *                 npem_get_rates calculates the rate matrix ("DM" in
 *                                the notation of Meng and Rubin)
 *                 npem_get_var   calculates the est'd var-cov matrix
 *
 *   npem_start.c: npem_start      program to get starting points
 *                 piksrt          simple sorting program from N.R. in C
 *                 get_seed        get seed for random number generator
 *                 gasdev          generate normal deviates
 *                 ran1            generate uniform(0,1) deviates
 *
 *   npem_matrix.c:ludcmp          LU decomposition
 *                 lubksb          forward-backward substitution
 *                 invert_matrix   matrix inversion
 *
 *********************************************************************/

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
 *                  or cells/well / 10^6 cells.)
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
 *                  (same length as n)
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
            int *n, int *ns, int *len_plate, int *len_groups,
            double *ests, double *k, double *ksq, double *loglik,
            int *maxit, double *tol, int *maxk, int *prnt,
            int *intwork, double *doublework,
            int *use_order_constraint);

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
              int *maxk);

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
              double *doublework, int *use_order_constraint);

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
            double *ests, double *loglik, int *maxk);

/*********************************************************************
 *
 * gammaln
 *
 * This function calculates the log of the gamma function. (Recall
 * that gamma(x+1) = x!.)
 *
 * Input: a double
 *
 * Output: a double
 *
 *********************************************************************/

double gammaln(double x);


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
             int *maxit);

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
                  double *ksq, double *infor, double *doublework);

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
 *      use_order_constraint
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
                   int *use_order_constraint, int *maxit);

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
 *     intwork = vector of length 2 * n_param to use as workspace
 *
 *     doublework = vector of length *** to use as workspace
 *
 **********************************************************************/

void npem_get_var(double *infor, double *rates, int *n_param,
                 double *var, double *se, int *intwork,
                 double *doublework);



void npem_start(double *y, int *n_plates, double *c, int *n_groups,
               int *n, int *ns, int *len_plate, double *cv,
               double *n_sd, double *ests, double *work);

void piksrt(int n, double *arr);

int get_seed(void);
double gasdev(int *idum);
double ran1(int *idum);



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
            int *errorflag);

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

void lubksb(double *a, int *n, int *indx, double *b);



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
                   double *doublework, int *errorflag);
