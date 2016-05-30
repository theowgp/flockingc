#include "st_descent.h"





int st_descent(
    double      grad_tol , /* StopRule = 1: |g|_infty <= max (grad_tol,
                                            StopFac*initial |g|_infty) [default]
                              StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
    double            *x , /* input: starting guess, output: the solution */
    int              dim , /* problem dimension (also denoted n) */
    double    (*cg_value)  /* user provided routine to return the function */
               (double *), /* value cg_value(x) at x */
    void       (*cg_grad)  /* user provided routine cg_grad (g, x), g is*/
     (double *, double *), /* the gradient at x*/
    double           step, /* initial step for line search */
    double              h
                        
)
{

  step = 1;

  double g[dim];
  double  s;  

  cg_grad(g, x); 
 

  int iter = 1;
  while(normL2(g, dim, h) > grad_tol && iter < 10000)
  {
    //s = determine_step_size();
    s = step;

    // update x: x +=  -s*g; 
    updatex(x, x, g, s, dim);

    // update gradient
    cg_grad(g, x); 

    printf("||g|| = %f, value = %f, iter = %d\n", norm(g, dim), cg_value(x), iter);

    iter++;

  }


  return 0;
}




void updatex(double *temp, double *x, double *g, double s, int n)
{  

  int i;
  for(i=0; i<n; i++)
  {
    temp[i] = x[i] - s*g[i];
  } 

 
}

