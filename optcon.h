#include "cg_descent.h"
typedef struct optcon_struct
{
    int n;    
    int nx;  
    int nc;  
    int ns; 
    double t_0; 
    double t_f; 
    double h;   
    double *state_0;
    double *state_f;
    double *state;
    double *costate;
    double *gx;
    double *gu;
    double *y;
    double *z;
    double *a;
    double *b;
    double *c;
    double *tau;
    double *cg_work;
} optcon_work;
typedef struct optcon_parameter
{
    int PrintLevel;
    int PrintFinal;
    int scheme;  
} optcon_parameter;
/* ---------------------------------------
    prototype of user defined functions 
   --------------------------------------- */
void (*optcon_f)(double *, double *, double *, double) ;// evaluate f(x)
void (*optcon_dphi)(double *, double *) ;// evaluate of dF/dx
double (*optcon_phi)(double *) ;//evaluate of F
void (*optcon_fu)(double *, double *, double *, double) ;//evaluate df/du
void (*optcon_fx)(double *, double *, double *, double) ;//evaluate df/dx
/* ---------------------------------------
       prototype of optcon functions 
   --------------------------------------- */
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
    double *a,    /*  Runge-Kutta scheme A   */
    double *b,    /*  Runge-Kutta scheme B   */
    optcon_parameter *parm
) ;
void optcon_gradient 
(
    double *g,      // ptr of the gradient vector of F respect to u
    double *control // ptr of the control
) ;
double optcon_cost
(
    double *control
) ;
void optcon_state 
( 
    double *control 
) ;
void optcon_costate 
( 
    double *control 
) ;
void optcon_combine 
( 
    double *g, 
    double *control 
) ;
int optcon
(
    double grad_tol,  /*  gradient tolerance to be used in cg_descent */
    int n,        /*  number of time mesh intervals*/
    int nx,       /*  number of state variables    */
    int nc,       /*  number of controllers        */
    int ns,       /*  number of stage in Runge-Kutta scheme*/
    double t_0,       /*     start time            */
    double t_f,       /*     finish time           */
    double *control,  /*  control variables  */
    double *state_0,   /*  initial states    */
    double *state,    /*  state variables    */
    double *a,    /*  Runge-Kutta scheme A   */
    double *b,     /*  Runge-Kutta scheme B   */
    double (*phi)(double *),
    void (*dphi)(double *, double *),
    void (*f)(double *, double *, double *, double),
    void (*fx)(double *, double *, double *, double),
    void (*fu)(double *, double *, double *, double)
) ;
