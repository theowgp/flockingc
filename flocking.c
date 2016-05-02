#include "flocking.h"
#include "params1.h"



//  f ////////////////////////////////////////////////////////////////////////////////////////////////////////

void S(double *res, double *x, int N, int d)
{
	res += (N+1)*d;
	x   += (N+1)*d;

	int i;
	for(i=0; i<(N+1); i++)
	{
		Si(res, x, d); 
		x  +=d;
		res+=d;
	}
}

void Si(double *res, double *v, int d)
{
	int i;

	for(i=0; i<d; i++)
	{
		res[i] = (alpha - beta*norm(v, d)*norm(v, d))*v[i];
	}
}




void fx(double *res, double *x, int N, int d)
{
	x += (N+1)*d;

	int i, j, k;
	k=0;
	for(i=0; i<=N; i++)
	{
		for (j=0; j<d; j++)
		{
			res[k] = x[k];
			k++;
		}
	}
}


void M(double *res, double *x, int N, int d)
{
	int i, j, k, l;

	for(i=0; i<(N+1); i++)
	{
		for(k=0; k<d; k++)
		{
			double temp1 = 0;
			for(j=0; j<(N+1); j++)
			{
				if(i!=j)
				{
					double temp2[2];
					for(l=0; l<d; l++)
					{	
						temp2[l] = x[i*d+l] - x[j*d+l];
					}
					
					double ro = norm(temp2, d);

					temp1 += temp2[k]*(Ca/la*exp(-ro/la) - Cr/lr*exp(-ro/lr))/ro;
				}
			}
			res[(N+1)*d + i*d + k] += -temp1/(N+1);
		}
	}
}

void L(double *res, double *x, int N, int d)
{
	int i, j, k, l;

	for(i=1; i<(N+1); i++)
	{
		for(k=0; k<d; k++)
		{
			double temp2[2];
			for(l=0; l<d; l++)
			{	
				temp2[l] = x[i*d+l] - x[l];
			}
					
			double ro = norm(temp2, d);

			res[(N+1)*d + i*d + k] += -gamma1*temp2[k]*(Ca0/la0*exp(-ro/la0) - Cr0/lr0*exp(-ro/lr0))/ro;
		}
	}
}


//  Gf ////////////////////////////////////////////////////////////////////////////////////////////////

void GS(double *res, double *x, int N, int d, int nx)
{
	int i, j, k, l;

	

	for(i=N+1; i<2*(N+1); i++)
	{
		double temp = norm(&x[i*d], d);

		for(k=i*d; k<i*d+d; k++)
		{
			for(l=i*d; l<i*d+d; l++)
			{
				if(k==l)
				{
					res[map(k, l, nx)] += (alpha - beta*temp*temp);
				}
				res[map(k, l, nx)] += -2*beta*x[k]*x[l];
			}	
		}
	}


	

}




 //  Destination /////////////////////////////////////////////////////////////////////////////////////////////////////////
double* xdes(double t, int d)
{
	double *res = (double*) malloc(d*sizeof(double));
	// double res[d];

	int k;
	for(k=0; k<d; k++)
	{
		res[k] = t;
	}
	return res;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
double norm(double *v, int n)
{
	double temp = 0;

	int i;
	for(i = 0; i<n; i++)
	{
		temp += v[i]*v[i];
	} 

	return sqrt(temp);
}

int map(int i, int j, int nx)
{
	return i*nx + j;
}
int* imap(int k, int nx)
{
	int *res =	(int*) malloc(2*sizeof(int));
	res[0] = k%nx;
	res[1] = k/nx;
	return res;
}
