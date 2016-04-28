
#include "flocking.h"
#include "optcon.h"
#include "params2.h"



void pad(double *f, int n);

void my_f(double *f, double *x, double *u, double time);// evaluate of f(x)
void my_dphi(double *dphi, double *x_f); // evaluate of dF/dx
double my_phi(double *x_f); //evaluate of F
void my_fu(double *fu, double *x, double *u, double time); //evaluate of df/du
void my_fx(double *fx, double *x, double *u, double time); //evaluate of df/dx

int N = 1;
int d = 2;
int n = 100;
// int nx = 2*(N+1)*d + 1; 
int nx = 9;
int nc = 2;
int ns = 3;







int main(void) 
{
	 // double x[] = {2, 2, 2, 2, 2, 2, 3, 1, 5, 6, 7, 8, 999 };
     // double x[] = {2, 2, 2, 2, 2, 2, 3, 1, 0, 5, 6, 0, 999 };
     double x[] = {2, 2, 2, 2, 3, 1, 5, 6, 999 };
     double u[] = {1, 1};
     // double f[nx];
     double *f;
     f = (double*) malloc(nx*sizeof(double));
     // f[nx-1] = 999;

     my_f(f, x, u, 0);

     pad(f, nx);

     // int testtemp1 = 7%2;
     // int testtemp2 = 7/2;
     // printf("\n k = %d\n", map(0, 2, 3));
     // int *testtemp = imap(6, 3);
     // printf(" i = %d\n j = %d\n\n", testtemp[0], testtemp[1]);
}



void pad(double *f, int n)
{
     printf("pda = \n");
     int i;
     for(i=0; i<nx; i++)
     {
       printf("%f  ", f[i]);
       // printf("%f  ", *f++);
     } 
     printf("\n");
}

void my_f(double *f, double *x, double *u, double time)// evaluate of f(x)
{
        
     // double testtemp;
	 // testtemp = norm(u, nc);

	 fx(f, x+(N+1)*d, N, d);          //What exactly is the vector x ??? Is it n*ns*nx+nx or just nx ???
	     
	 S(f+(N+1)*d, x+(N+1)*d, N, d);

	 int i;
	 for(i=0; i<d; i++)
	 {
	 	*(f + (N+1)*d + i) += u[i];
	 }
	    
	 *(f+nx-1) = 0.5*nu*norm(u, nc);  //What exactly is the vector u ??? Is it n*ns*nc    or just nc ???

	    
	 return;
 }


void my_dphi(double *dphi, double *x_f) // evaluate of dF/dx
{
     dphi[0] = 0;
     dphi[1] = .5;
     return;
     ///* delete after compile */
     //double a;
     //a = x_f[0];
}

double my_phi(double *x_f) //evaluate of F
{
     return .5*x_f[1];
}

void my_fu(double *fu, double *x, double *u, double time) //evaluate of df/du
{
     fu[0] = 1;
     fu[1] = 2*u[0];
     return;
     ///*delete after compile*/
     //double a, b;
     //a = x[0];
     //b = time;
}   

void my_fx(double *fx, double *x, double *u, double time) //evaluate of df/dx
{
     fx[0] = 0.5;
     fx[2] = 4*x[0];
     fx[1] = 0;
     fx[3] = 0;

     return;
}     
