
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

int N = 9;
int d = 2;
int n = 10;
// int nx = 2*(N+1)*d + 1; 
int nx = 41;
int nc = 2;
int ns = 3;

int nMPC = 1;



int main(void) 
{
	double t0, tf;

	double dt;

	t0 = 0;
	tf = T;
	dt = tf/nMPC;

	double state0[nx];
	double grad_tol = 1.0e-7;

	double *a, *b;
	int i, j;
	
	// double  states    [nMPC] [n*ns*nx+nx];
	// double  controls  [nMPC] [n*ns*nc];

	double  *states    [nMPC];
	double  *controls  [nMPC];
	
	// double  **states;
	// double  **controls;
	// states = (double**) malloc(nMPC*sizeof(double*));
	// controls = (double**) malloc(nMPC*sizeof(double*));
   
	for(i=0; i<nMPC; i++)
	{          
		states[i] =   (double*) malloc((n*ns*nx+nx)*sizeof(double));
		controls[i] = (double*) malloc(n*ns*nc*sizeof(double));
	}

	
	
	
	
	

	/*initial control*/
	for(i = 0; i < nMPC; i ++)
	{
		for ( j = 0; j < n*ns*nc; j ++ ) 
		{
			controls[i][j] = 0;
		}

	}
	
	/*initial state0*/
	// state0[0] = 0;
	// state0[1] = 0;
	// state0[2] = 1;
	// state0[3] = 0;
	// state0[4] = 0;
	// state0[5] = 0.5;
	// state0[6] = 0;
	// state0[7] = 0.5;
	// state0[8] = 0;

	// initialize state
	state0[0] = -1;
	state0[1] = 0;
	state0[(N+1)*d + 0] = 0;
	state0[(N+1)*d + 1] = -1;

	state0[d + 0] = 0;
	state0[d + 1] = 0;
	state0[(N+1)*d + d + 0] = 0;
	state0[(N+1)*d + d + 1] = -0.2;
	
	int k;
	for (i=2; i<N+1; i++)
	{
		for(k=0; k<d; k++)
		{
			state0[i*d + 0] = state0[(i-1)*d + 0] + 0.01;     
			state0[i*d + 1] = 0;
		}
		for(k=0; k<d; k++)
		{
			state0[(N+1)*d + i*d + 0] =  0;     
			state0[(N+1)*d + i*d + 1] =  -0.2;
		}
	}
	 

	double t1, t2;
	t1 = t0;
	t2 = dt;
	T  = t2;
	// optcon(grad_tol, n, nx, nc, ns, t1, t2, &controls[0][0], &state0, 
	//        &states[0][0], a, b, my_phi, my_dphi, my_f, my_fx, my_fu);
	// for(i = 1; i < nMPC; i ++)
	// {
	//      t1 += dt;
	//      t2 += dt;
	//      optcon(grad_tol, n, nx, nc, ns, t1, t2, &controls[i][0], &states[i-1][n*ns*nx], 
	//        &states[i][0], a, b, my_phi, my_dphi, my_f, my_fx, my_fu);

	// }
	optcon(grad_tol, n, nx, nc, ns, t1, t2, controls[0], &state0, 
		  states[0], a, b, my_phi, my_dphi, my_f, my_fx, my_fu);
	for(i = 1; i < nMPC; i ++)
	{
		t1 += dt;
		t2 += dt;
		T  =  t2; 
		optcon(grad_tol, n, nx, nc, ns, t1, t2, controls[i], (states[i-1] + n*ns*nx), 
		  states[i], a, b, my_phi, my_dphi, my_f, my_fx, my_fu);

	}
	
	
	
	FILE *fs = fopen("state.txt", "w");
    wtf(states, fs, n, ns, nx, N, d, nMPC, t0, tf);    
	fclose(fs);

	FILE *fc = fopen("control.txt", "w");
    wtf(controls, fc, n, ns, nc, N, d, nMPC, t0, tf);    
	fclose(fc);
    
	
	printf("\n");









  //     // double x[] = {2, 2, 2, 2, 2, 2, 3, 1, 5, 6, 7, 8, 999 };
  //    // double x[] = {2, 2, 2, 2, 2, 2, 3, 1, 0, 5, 6, 0, 999 };
  //    // double x[] = {1, 0, 2, 0, 1, 2, 3, 4, 999 };
	 // double *x;
	 // x = state0;
	 // // double x[] = {1, 1, 1, 2, 999 };
  //    double u[] = {1, 2};
	


	// // double *f;
	// // f = (double*) malloc(nx*sizeof(double));  // what is the difference ???????????
	// double f[nx];                                // what is the difference ???????????
	// my_f(f, x, u, 0);
	// printf("\nf = \n");
	// pad(f, nx);

	// // int testtemp1 = 7%2;
	// // int testtemp2 = 7/2;
	// // printf("\n k = %d\n", map(0, 2, 3));
	// // int *testtemp = imap(6, 3);
	// // printf(" i = %d\n j = %d\n\n", testtemp[0], testtemp[1]);

	// // double *matrixf;
	// // matrixf = (double*) malloc(nx*nx*sizeof(double));
	// double matrixf[nx*nx];
	// my_fx(matrixf, x, u, 0);
	// printf("\nfx =");
	// pmd(matrixf,  nx, nx);
	// // pmd(matrixf, nx*nx, 1);
    

	// // double *matrixu;
	// // matrixu = (double*) malloc(nx*d*sizeof(double));
	// double matrixu[nx*d];
	// my_fu(matrixu, x, u, 0);
	// printf("\nfu =");
	// pmd(matrixu,  nx, d);
	// // pmd(matrixu, nx*d, 1);
    

	// double valuephi;
	// valuephi = my_phi(x);
	// printf("\nphi = %f\n", valuephi);

	// // double *vectordphi;
	// // vectordphi = (double*) malloc(nx*sizeof(double));
	// double vectordphi[nx];
	// my_dphi(vectordphi, x);
	// printf("\ndphi =");
	// pmd(vectordphi,  nx, 1);
    

	// double valuel1 = l1(x, N, d);
	// printf("\nl1 = %f\n", valuel1);



	// printf("\n");
}







void my_f(double *f, double *x, double *u, double time)// evaluate of f(x)
{
	   
	// double testtemp;
	 // testtemp = norm(u, nc);

	 fx(f, x, N, d);          //What exactly is the vector x ??? Is it n*ns*nx+nx or just nx ???
		
	 S(f, x, N, d);

	 M(f, x, N, d);

	 L(f, x, N, d);

	 int i;
	 for(i=0; i<d; i++)
	 {
		*(f + (N+1)*d + i) += u[i];
	 }
	    
	 *(f+nx-1) = 0.5*nu*norm(u, nc)*norm(u, nc) + l1(x, N, d);  //What exactly is the vector u ??? Is it n*ns*nc    or just nc ???

	    
	 return;
}


void my_dphi(double *dphi, double *x_f) // evaluate of dF/dx
{
	 setmt0(dphi, nx, d);

	 
	 int k;
	 for(k=0; k<d; k++)
	 {
		dphi[k] = x_f[k] - xdes(T, k);
	 }
	 

	 dphi[nx-1] = 1;

	return;
}

double my_phi(double *x_f) //evaluate of F
{
	 double temp[d];
	 

	 int k;
	 for(k=0; k<d; k++)
	 {
		temp[k] = x_f[k] - xdes(T, k);
	 }
	 

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

	 GM(fx, x, N, d, nx);
	 
	 GL(fx, x, N, d, nx);

	 Gl1(fx, x, N, d, nx);

	return;
}     












