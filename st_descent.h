#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#include "space.h"

int st_descent /*  return  0 (convergence tolerance satisfied)
                           1 (change in func <= feps*|f|)
                           2 (total iterations exceeded maxit)
                           3 (slope always negative in line search)
                           4 (number secant iterations exceed nsecant)
                           5 (search direction not a descent direction)
                           6 (line search fails in initial interval)
                           7 (line search fails during bisection)
                           8 (line search fails during interval update)
                           9 (debugger is on and the function value increases)
                          -1 (parameter file not found)
                          -2 (missing parameter value in parameter file)
                          -3 (comment in parameter file too long) */
(
    double      grad_tol , /* StopRule = 1: |g|_infty <= max (grad_tol,
                                            StopFac*initial |g|_infty) [default]
                              StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
    double            *x , /* input: starting guess, output: the solution */
    int              dim , /* problem dimension (also denoted n) */
    double    (*cg_value)  /* user provided routine to return the function */
               (double *), /* value cg_value(x) at x */
    void       (*cg_grad)  /* user provided routine cg_grad (g, x), g is*/
     (double *, double *), /* the gradient at x*/
    // double         *work , /* working array with at least 4n elements */
    double          step,   /* initial step for line search */
    double             h
                        
);



void updatex(double *temp, double *x, double *g, double s, int n);