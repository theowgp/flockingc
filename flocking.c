#include "flocking.h"
#include "params1.h"





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