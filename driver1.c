
#include "flocking.h"
#include "optcon.h"
#include "params2.h"
#include "plugin1.h"




void my_f(double *f, double *x, double *u, double time);// evaluate of f(x)
void my_dphi(double *dphi, double *x_f); // evaluate of dF/dx
double my_phi(double *x_f); //evaluate of F
void my_fu(double *fu, double *x, double *u, double time); //evaluate of df/du
void my_fx(double *fx, double *x, double *u, double time); //evaluate of df/dx





double T = 1;

int N = 0;
int d = 2;
int n = 100;
// int nx = 2*(N+1)*d + 1; 
int nx = 5;
int nc = 2;
int ns = 3;






int main(void) 
{
	 double t0, tf ;
     double state0[nx];
     double grad_tol = 1.0e-7;
     double *state;
     double *control;
     double *a, *b;
     int i;
     
     t0 = 0;
     tf = T;
     /* allocate memory */
     state = (double*) malloc((n*ns*nx+nx)*sizeof(double));
     control = (double*) malloc(n*ns*nc*sizeof(double));
     /*initial control*/
     for ( i = 0; i < n*ns*nc; i ++ ) 
     {
         // control[i] = .5;
     	control[i] = 0;
     }
     /*initial state0*/
     state0[0] = 0;
     state0[1] = 0;
     state0[2] = 0;
     state0[3] = 0.5;
     state0[4] = 0;
     /*call optcon_xrk*/
     optcon(grad_tol, n, nx, nc, ns, t0, tf, control, &state0[0], 
            state, a, b, my_phi, my_dphi, my_f, my_fx, my_fu);
     
     free(state);
     free(control);

     printf("\n");
}







void my_f(double *f, double *x, double *u, double time)// evaluate of f(x)
{
        
     // double testtemp;
	 // testtemp = norm(u, nc);

	 fx(f, x, N, d);          //What exactly is the vector x ??? Is it n*ns*nx+nx or just nx ???
	     
	 S(f, x, N, d);

	 int i;
	 for(i=0; i<d; i++)
	 {
	 	*(f + (N+1)*d + i) += u[i];
	 }
	    
	 *(f+nx-1) = 0.5*nu*norm(u, nc)*norm(u, nc);  //What exactly is the vector u ??? Is it n*ns*nc    or just nc ???

	    
	 return;
}


void my_dphi(double *dphi, double *x_f) // evaluate of dF/dx
{
 	 setmt0(dphi, nx, d);

 	 double *des = xdes(T, d);

 	 int k;
	 for(k=0; k<d; k++)
	 {
	 	dphi[k] = x_f[k] - des[k];
	 }
	 free(des);

	 dphi[nx-1] = 1;

     return;
}

double my_phi(double *x_f) //evaluate of F
{
	 double *temp =	(double*) malloc(d*sizeof(double));
	 double *des = xdes(T, d);

	 int k;
	 for(k=0; k<d; k++)
	 {
	 	temp[k] = x_f[k] - des[k];
	 }
	 free(des);

     return 0.5*norm(temp, d)*norm(temp, d) + x_f[nx-1];
}

void my_fu(double *fu, double *x, double *u, double time) //evaluate of df/du
{ 
     setmt0(fu, nx, d);

     int i, k, l;

     i = (N+1);
     for(k=0; k<d; k++)
	 {
		 fu[map(i*d + k, k, d)] = 1;
	 }

	 i = 2*(N+1);
	 for(k=0; k<d; k++)
	 {
		 fu[map(i*d, k, d)] = nu*u[k];
	 }
	 

     
     return;
}   

void my_fx(double *fx, double *x, double *u, double time) //evaluate of df/dx
{

	   
	 setmt0(fx, nx, nx);


	 // set the dfx/dv part which is an identity matrix
	 int i, k, l;
	 for(i=0; i<N+1; i++)
	 {
		 for(k=i*d; k<i*d+d; k++)
		 {
			 fx[map(k, k + (N+1)*d, nx)] = 1;
						
		 }
	 }

	 GS(fx, x, N, d, nx);

     return;
}     












