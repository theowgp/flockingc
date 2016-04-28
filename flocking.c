#include "flocking.h"
#include "params1.h"



//  f ////////////////////////////////////////////////////////////////////////////////////////////////////////

void S(double *res, double *v, int N, int d)
{
	int i;

	for(i=0; i<=N; i++)
	{
		Si(res, v, d); 
		v  +=d;
		res+=d;
	}
}

void Si(double *res, double *v, int d)
{
	int i;

	for(i=0; i<d; i++)
	{
		res[i] =   (alpha - beta*norm(v, d))*v[i];
	}
}




void fx(double *res, double *v, int N, int d)
{
	int i, j, k;

	k=0;
	for(i=0; i<=N; i++)
	{
		for (j=0; j<d; j++)
		{
			res[k] = v[k];
			k++;
		}
	}
}



//  Gf ////////////////////////////////////////////////////////////////////////////////////////////////

void GS(double *res, double *v, int N, int d)
{
	int i;

	for(i=0; i<=N; i++)
	{
		Si(res, v, d); 
		v  +=d;
		res+=d;
	}
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
	return j*nx + i;
}
int* imap(int k, int nx)
{
	int *res =	(int*) malloc(2*sizeof(int));
	res[0] = k%nx;
	res[1] = k/nx;
	return res;
}
