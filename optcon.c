#include "optcon.h"
/* 
       ________________________________________________________________
      |               C-code Version 1.0  (May 5, 2006)                |
      |                William W. Hager  and   Shuo Li                 |
      |               hager@math.ufl.edu    shuo76@ufl.edu             |
      |                   Department of Mathematics                    |
      |                     University of Florida                      |
      |                 Gainesville, Florida 32611 USA                 |
      |                      352-392-0281 x 244                        |
      |                                                                |
      |                     Copyright by Shuo Li                       |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/Control         |
      |________________________________________________________________|

   optcon is a routine that solves an unconstrained control problem:

            min phi(x(t_f))

            subject to x' = f(x, u, t), x(t_0) = x_0

   The continuous problem is discretized by an explicit Runge-Kutta scheme:

       x_{k+1} = x_k + h \sum_{i=1}^s b_i f(y_i, u_{ki})

       y_i = x_k + h \sum_{j=1}^s a_{ij} f(y_j, u_{kj})

   The optimization is done using the conjugate gradient scheme CG_DESCENT.
       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|
       _________________________________________________________________
      |Note: The file optcon_c.parm must be placed in the directory    |
      |      where the code is run                                     |
      |________________________________________________________________|

   In order to calculate the discretized gradient, the user should
   define and pass the following functions to optcon:

   1. void my_f(double *f, double *x, double *u, double time)
           *x, *u, and time are inputs, *f is output.
   2. void my_dphi(double *dphi, double *x_f)
           *x_f is input, *dphi is output.
   3. double my_phi(double *x_f)
           *x_f is input, the output is the evaluation of phi, 
           which has double precision.
   4. void my_fu(double *fu, double *x, double *u, double time)
           *x, *u, and time are inputs, *fu is output.
   5. void my_fx(double *fx, double *x, double *u, double time)
           *x, *u, and time are inputs, *fx is output.

   *************************************************************
   The inputs and output of optcon are:
   Inputs:
         grad_tol - gradient tolerance to be used in cg_descent
         n        - number of time mesh intervals
         nx       - number of state variables
         nc       - number of controllers
         ns       - number of stage(s) in Runge-Kutta scheme
         t_0      - initial time
         t_f      - final time
         *control - pointer points to the initial control
         *state_0 - pointer points to the initial state
         *state   - pointer points to the state
         *a       - pointer points to the R-K matrix A
         *b       - pointer points to the R-K vector b
         double (*phi)(double *) - pointer points to the user routine 
                                   that evaluates terminal cost phi.
         void (*dphi)(double *, double *) - pointer points to the user
                      routine that evaluates gradient of the terminal cost.
         void (*f)(double *, double *, double *, double) - pointer points 
                                 to the user routine evaluates f.
         void (*fx)(double *, double *, double *, double) - pointer points 
                           to the user routine evaluates df/dx.
         void (*fu)(double *, double *, double *, double) - pointer points 
                           to the user routine evaluates df/du.
   Output:
          optcon returns an integer has the following meanning:
           0    - convergence tolerance satisfied
           1    - change in func <= feps*|f|
           2    - total CG iterations exceeded maxit
           3    - slope always negative in line search
           4    - number secant iterations exceed nsecant
           5    - search direction not a descent direction
           6    - line search fails in initial interval
           7    - line search fails during bisection
           8    - line search fails during interval update
           9    - debugger is on and the function value increases
          -1    - CG parameter file not found
          -2    - missing parameter value in CG parameter file
          -3    - comment in CG parameter file too long
          -4    - not enough memory available
          -5    - optcon_c.parm is missing
          -6    - 'PrintLevel' parameter is missing in optcon_c.parm
          -7    - 'PrintFinal' parameter is missing in optcon_c.parm
          -8    - 'scheme' parameter is missing in optcon_c.parm
          -9    - comment is too long in optcon_c.parm
          -10   - some b[i] = 0
          -11   - given implicit RK schemes
          -12   - scheme parameter is out of range
          -13   - t_f is less than t_0
 */
optcon_work  work;
int optcon
(
    double grad_tol,
    int n,
    int nx,
    int nc,
    int ns,
    double t_0,
    double t_f,
    double *control,
    double *state_0,
    double *state,
    double *a,
    double *b,
    double (*phi)(double *),
    void (*dphi)(double *, double *),
    void (*f)(double *, double *, double *, double),
    void (*fx)(double *, double *, double *, double),
    void (*fu)(double *, double *, double *, double)
)
{
    extern optcon_work work ;
    cg_stats Stats;
    optcon_parameter parm;
    double step;
    int setup_status, cg_status;
    int i;

    
    
    /* set the user-provide functions as extern functions */
    optcon_phi  = phi ;
    optcon_dphi = dphi ;
    optcon_f    = f ;
    optcon_fx   = fx ;
    optcon_fu   = fu ;
    /* setup the problem: 
         1. malloc the working space;
         2. initiate the global variable optcon_work;
         3. read parameters from optcon_c.parm;
         4. check the explicity of given RK scheme;
    */
    setup_status = optcon_setup( n, nx, nc, ns, t_0, t_f, 
                 state_0, state, a, b, &parm );
    switch (setup_status) {
        case -4:
        {
            printf ("No enough memory available. \n");
            return (-4) ;
        }
        case -5:
        {
            printf ("optcon parameter file (optcon_c.parm) "
                   "not found. \n") ;
            return (-5) ;
        }
        case -6:
        {
            printf ("'PrintLevel' parameter is missing in optcon_c.parm. \n") ;
            return (-6) ;
        }
        case -7:
        {
            printf ("'PrintFinal' parameter is missing in optcon_c.parm. \n") ;
            return (-7) ;
        }
        case -8:
        {
            printf ("'scheme' parameter is missing in optcon_c.parm. \n") ;
            return (-8) ;
        }       
        case -9:
        {
            printf ("the length of a comment statement in the "
                   "file cg_descent_c.parm exceeds the maximum allowed "
                   "length %i . \n", MAXLINE) ;
            return (-9) ;
        }
        case -10:
        {
            printf("some b[i] vanishes. This program only "
                  "handles the case of all b[i] nonzero. \n");
            return (-10) ;
        }
        case -11:
        {
            printf("This program treats explicit Runge-Kutta "
                  "schemes only. The given scheme is implicit. \n");
            return (-11) ;
        }
        case -12:
        {
            printf ("'scheme' parameter is out of range.\n") ;
            return (-12) ;
        }
        case -13:
        {
            printf ("final time should be greater than initial time. \n");
            return (-13) ;
        }
    }
    
    // cg_status = cg_descent (grad_tol, control, work.n*work.ns*work.nc, 
    //           optcon_cost, optcon_gradient, work.cg_work, step, &Stats);
    // if ( parm.PrintFinal || parm.PrintLevel ) 
    // {
    //      printf ("Final state: \n");
    //      for ( i = 0; i < nx; i ++ ) {
    //          printf ("final[%d] = %.10e\n", i+1, work.state_f[i]);
    //      }
    //      printf ("Terminal cost: %.10e\n", Stats.f);
    // }
    
    double h = (t_f-t_0)/n/ns;
    cg_status = st_descent (grad_tol, control, work.n*work.ns*work.nc, 
              optcon_cost, optcon_gradient, step,  h) ;

    free (work.cg_work);
    return (cg_status);
}
/*
    Function optcon_setup()
     input:
           n,             number of time mesh intervals
           nx,            number of state variables
           nc,            number of controllers
           ns,            number of stage in Runge-Kutta scheme
           t_0,           start time
           t_f,           finish time
           *state_0       ptr points to the initial state
           *state         ptr points to the state
           *a,            a pointer points to the Runge-Kutta a
           *b,            a pointer points to the Runge-Kutta b
           *parm          a pointer points to the optcon_parameter
    In parm file optcon_c.parm has the following
    PARAMETERS:
 
    PrintLevel - print the result from each iteration
    PrintFinal - print the final solution
    scheme - choose the Runge-Kutta scheme you want to use:
    0:          user provides
    
    1:               | 0  0 |            |1/2|
                A =  | 1  0 |        b = |1/2|
    
    2:              | 0  0  0|           |1/6|
                A = |.5  0  0|       b = |2/3|
                    |-1  2  0|           |1/6|
    
    3:              | 0  0  0|           |2/9|
                A = |.5  0  0|       b = |1/3|
                    |0 .75  0|           |4/9|
    
    4:              | 0  0  0|           |1/3|
                A = |.5  0  0|       b = |1/3|
                    |.5 .5  0|           |1/3|
    
    5:              | 0   0   0|         |1/6|
                A = | 1   0   0|     b = |1/6|
                    |.25 .25  0|         |2/3|
    
    6:              | 0  0  0  0|        |1/6|
                A = |.5  0  0  0|    b = |1/3|
                    | 0 .5  0  0|        |1/3|
                    | 0  0  1  0|        |1/6|
    DEFAULT PARAMETER VALUES:
    PrintLevel: 0
    PrintFianl: 1
    scheme: 2
    output: an integer from 0 to -13:
            
             0    problem initialized successfully
            -4    no enough memory available
            -5    optcon_c.parm is missing
            -6    'PrintLevel' parameter is missing in optcon_c.parm
            -7    'PrintFinal' parameter is missing in optcon_c.parm
            -8    'scheme' parameter is missing in optcon_c.parm
            -9    comment is too long
            -10   some b[i] = 0
            -11   given implicit RK schemes
            -12   scheme parameter is out of range
            -13   t_f is less than t_0
    The function will do:
        1: initialize structure optcon_work
        2: initialize structure optcon_parameter
        3: the memory will be evaluated and allocated for the problem
        4: the given RK scheme will be checked for explicity
*/
int optcon_setup 
(
    int n,        /*  number of time mesh intervals*/
    int nx,       /*  number of state variables    */
    int nc,       /*  number of controllers        */
    int ns,       /*  number of stage in Runge-Kutta scheme*/
    double t_0,   /*     start time            */
    double t_f,   /*     finish time           */
    double *state_0,  /*  initial states    */
    double *state,    /*  state variables    */
    double *a,    /*  Runge-Kutta matrix A   */
    double *b,    /*  Runge-Kutta vector b   */
    optcon_parameter *parm
)
{
    extern optcon_work  work ;
    FILE *ParmFile ; 
    char junk [MAXLINE+1] ; 
    double t;         
    int Maxsize;  
    int info;
    int i, j;
    /* define the working space */
    Maxsize = 4*n*ns*nc + nx*nx + nx*nc + 2*ns*ns + 2*ns + 2*nx + n*ns*nx;
    /* malloc computer memory */
    work.cg_work = (double*) malloc(Maxsize*sizeof(double));
    /* if there is not enough memory, return error */
    if ( work.cg_work == NULL ) return (-4) ;
    /* work.cg_work(used in CG_DESCENT) needs 4*n*ns*nc elements*/
    work.tau = work.cg_work + 4*n*ns*nc;/* work.tau needs ns elements */
    work.c = work.tau + ns;/* work.c needs ns*ns elements */
    work.b = work.c + ns*ns;/* work.b needs ns elements */
    work.a = work.b + ns;/* work.a needs ns*ns elements */
    work.z = work.a + ns*ns;/* work.z needs nx elements */
    work.y = work.z + nx;/* work.y needs nx elements */
    work.gu = work.y + nx;/* work.gu needs nx*nc elements */
    work.gx = work.gu + nx*nc;/* work.gx needs nx*nx elements */
    work.costate = work.gx + nx*nx;/* work.costate needs n*ns*nx elements */
    work.state = state;
    work.state_f = state + n*ns*nx;
    work.state_0 = state_0;
    work.ns = ns;
    work.nc = nc;
    work.nx = nx;
    work.n = n;
    if ( t_f <= t_0 ) return (-13) ;
    work.t_0 = t_0;
    work.t_f = t_f;
    work.h = ( t_f - t_0 )/n ;
    ParmFile = fopen ("optcon_c.parm", "r") ;
    if ( ParmFile == NULL ) return (-5) ;
    info = fscanf ( ParmFile, "%i", &( parm->PrintLevel ) ) ;
    if ( info != 1 )
    {
        return (-6) ;
    }
    fgets (junk, MAXLINE, ParmFile) ;
    if (strlen (junk) >= MAXLINE-1) return (-9) ;
    info = fscanf ( ParmFile, "%i", &( parm->PrintFinal ) ) ;
    if ( info != 1 )
    {
        return (-7) ;
    }
    fgets (junk, MAXLINE, ParmFile) ;
    if (strlen (junk) >= MAXLINE-1) return (-9) ;
    info = fscanf ( ParmFile, "%i", &( parm->scheme ) ) ;
    if ( info != 1 )
    {
        return (-8) ;
    }
    fgets (junk, MAXLINE, ParmFile) ;
    if (strlen (junk) >= MAXLINE-1) return (-9) ;
    fclose (ParmFile) ;
    /*  initialize Runge-Kutta scheme 
        if user chose to use his/her own scheme, 
        the explicity needs to be checked
    */
    if (parm->scheme == 1) {
                work.a[0] = 0.;
                work.a[1] = 0.;
                work.a[2] = 1.;
                work.a[3] = 0.;
                work.b[0] = .5;
                work.b[1] = .5;
    }
    else if (parm->scheme == 2) {
                work.a[0] = 0.;
                work.a[1] = 0.;
                work.a[2] = 0.;
                work.a[3] = .5;
                work.a[4] = 0.;
                work.a[5] = 0.;
                work.a[6] = -1.;
                work.a[7] = 2.;
                work.a[8] = 0.;
                work.b[0] = 1./6.;
                work.b[1] = 2./3.;
                work.b[2] = 1./6.;
    }
    else if (parm->scheme == 3) {
                work.a[0] = 0.;
                work.a[1] = 0.;
                work.a[2] = 0.;
                work.a[3] = .5;
                work.a[4] = 0.;
                work.a[5] = 0.;
                work.a[6] = 0.;
                work.a[7] = .75;
                work.a[8] = 0.;
                work.b[0] = 2./9.;
                work.b[1] = 1./3.;
                work.b[2] = 4./9.;
    }
    else if (parm->scheme == 4) {
                work.a[0] = 0.;
                work.a[1] = 0.;
                work.a[2] = 0.;
                work.a[3] = .5;
                work.a[4] = 0.;
                work.a[5] = 0.;
                work.a[6] = .5;
                work.a[7] = .5;
                work.a[8] = 0.;
                work.b[0] = 1./3.;
                work.b[1] = 1./3.;
                work.b[2] = 1./3.;
    }
    else if (parm->scheme == 5) {
                work.a[0] = 0.;
                work.a[1] = 0.;
                work.a[2] = 0.;
                work.a[3] = 1.;
                work.a[4] = 0.;
                work.a[5] = 0.;
                work.a[6] = .25;
                work.a[7] = .25;
                work.a[8] = 0.;
                work.b[0] = 1./6.;
                work.b[1] = 1./6.;
                work.b[2] = 2./3.;
    }
    else if (parm->scheme == 6) {
                work.a[0] = 0.;
                work.a[1] = 0.;
                work.a[2] = 0.;
                work.a[3] = 0.;
                work.a[4] = .5;
                work.a[5] = 0.;
                work.a[6] = 0.;
                work.a[7] = 0.;
                work.a[8] = 0.;
                work.a[9] = .5;
                work.a[10] = 0.;
                work.a[11] = 0.;
                work.a[12] = 0.;
                work.a[13] = 0.;
                work.a[14] = 1.;
                work.a[15] = 0.;
                work.b[0] = 1./6.;
                work.b[1] = 1./3.;
                work.b[2] = 1./3.;
                work.b[3] = 1./6.;
    }    
    else if (parm->scheme == 0)
    {
        /* -------------------------
          check the given schemes  
        ------------------------- */
        for ( i = 0; i < ns; i ++ )
        {
            if ( b[i] == 0 )
            {
               return(-10); /* b[i] vanished */
            }
            work.b[i] = b[i];
            for ( j = 0; j < ns; j ++ )
            {
                if ( i >= j )
                {
                   if ( a[ns*j+i] != 0 )
                   {
                        /* implicit schemes */
                        return (-11);
                   }
                }
                work.a[i*ns+j] = a[i*ns+j];
            }
        }
    }
    else return (-12); // "scheme" parameter out of range.
    /* initiate tau and c */
    for ( i = 0; i < ns; i ++ )
    {
        t = 0;
        for ( j = 0; j < ns; j ++ )
        {
            t = t + work.a[ns*i+j];
            /* initialize c */
            work.c[i*ns+j] = work.b[j]*work.a[ns*j+i]/work.b[i];
        }
        /* initialize tau */
        work.tau[i] = t;
    }
    return (0);
}
/*
   optcon_gradient - evaluates gradient of cost respect to u.
   input:
         *control- ptr points to the control array
   output:
         *g      - ptr points to the gradient array
*/
void optcon_gradient 
(
    double *g,      // ptr of the gradient vector of F respect to u
    double *control // ptr of the control
)
{
    optcon_state ( control );
    optcon_costate ( control );
    optcon_combine ( g, control );
    return ;
}
/*
   optcon_value - evaluates phi(x(t_f)).
   input:
         *control- ptr points to the control array
   output:
          return the value of cost
*/
double optcon_cost
(
    double *control
)
{
    extern optcon_work work;
    double f;
    optcon_state ( control );
    f = optcon_phi ( work.state_f );
    return f ;
}
void optcon_state ( double *control ) 
{
    extern optcon_work work;
    int i, j, k, l;
    int kx, kc, ix, jx, jc;
    double time;
    for ( i = 0; i < work.nx; i ++ )
    {
        work.z[i] = work.state_0[i];
    }
    for ( k = 0; k < work.n; k ++ )
    {
        time = work.t_0 + work.h*k;
        kx = k*work.ns*work.nx;
        kc = k*work.ns*work.nc;
        for ( i = 0; i < work.ns; i ++ )
        {
            ix = kx + i*work.nx;
            for ( l = 0; l < work.nx; l ++ )
            {
                work.state[ix+l] = work.z[l];
            }
        }
        for ( j = 0; j < work.ns; j ++ )
        {
            jx = kx + j*work.nx;
            jc = kc + j*work.nc;
            optcon_f(work.y, &work.state[jx], &control[jc], 
                     time+work.h*work.tau[j]);
            for ( i = j+1; i < work.ns; i ++ )
            {
                ix = kx + i*work.nx;
                for ( l = 0; l < work.nx; l ++ )
                {
                    work.state[ix+l] = work.state[ix+l] + 
                                     work.h*work.a[i*work.ns+j]*work.y[l];
                }
            }
            for ( l = 0; l < work.nx; l ++ )
            {
                work.z[l] = work.z[l] + work.h*work.b[j]*work.y[l];
            }
        }
    }
    for ( l = 0; l < work.nx; l ++ )
    {
        work.state_f[l] = work.z[l];
    }
    return;
}    
void optcon_costate ( double *control )
{
     extern optcon_work work;
     int i, j, k, l;
     int ix, kx, kc, jx, jc;
     double time;
     double s;
     /*call dphi update z in Work*/
     optcon_dphi ( work.z, work.state_f ) ;
     for ( k = work.n-1; k > -1; k -- )
     {
         kx = k*work.ns*work.nx;
         kc = k*work.ns*work.nc;
         for ( i = 0; i < work.ns; i ++ )
         {
             ix = kx + i*work.nx;
             for ( l = 0; l < work.nx; l ++ )
             {
                 work.costate[ix+l] = work.z[l];
             }
         }
         time = work.t_0 + work.h*k;
         for ( j = work.ns-1; j > -1; j -- )
         {
             jx = kx + j*work.nx;
             jc = kc + j*work.nc;
             optcon_fx(work.gx, &work.state[jx], &control[jc], 
                       time + work.h*work.tau[j]);
             for ( l = 0; l < work.nx; l ++ )
             {
                 s = 0;
                 for ( i = 0; i < work.nx; i ++ )
                 {
                     s = s + work.costate[jx+i]*work.gx[i*work.nx+l];
                 }
                 work.y[l] = s;
             }
             for ( i = 0; i < j; i ++ )
             {
                 ix = kx + i*work.nx;
                 for ( l = 0; l < work.nx; l ++ )
                 {
                     work.costate[ix+l] = work.costate[ix+l] + 
                                        work.h*work.c[i*work.ns+j]*work.y[l];
                 }
             }
             for ( l = 0; l < work.nx; l ++ )
             {
                 work.z[l] = work.z[l] + work.h*work.b[j]*work.y[l];
             }
         }
     }
     return ;
}
void optcon_combine ( double *g, double *control )
{
     extern optcon_work work;
     int i, j, k, l;
     int kx, kc, jx, jc;
     double time;
     double s;
     for ( k = 0; k < work.n; k ++ )
     {
         kx = k*work.ns*work.nx;
         kc = k*work.ns*work.nc;
         time = work.t_0 + work.h*k ;
         for ( j = 0; j < work.ns; j ++ )
         {
             jx = kx + j*work.nx;
             jc = kc + j*work.nc;
             optcon_fu(work.gu, &work.state[jx], &control[jc], 
                       time + work.h*work.tau[j]);
             for ( i = 0; i < work.nc; i ++ )
             {
                 s = 0;
                 for ( l = 0; l < work.nx; l ++ )
                 {
                     s = s + work.costate[jx+l]*work.gu[l*work.nc+i];
                 }
                 g[jc+i] = work.h*work.b[j]*s;
             }
         }
     }
     return ;
}
