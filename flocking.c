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
					double temp2[d];
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
			double temp2[d];
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
	

	for(i=0; i<N+1; i++)
	{
		double temp = norm(&x[(N+1)*d + i*d], d);

		for(k=0; k<d; k++)
		{
			for(l=0; l<d; l++)
			{
				if(k==l)
				{
					res[map((N+1)*d + i*d + k, (N+1)*d + i*d + l, nx)] += (alpha - beta*temp*temp);
				}
				res[map((N+1)*d + i*d + k, (N+1)*d + i*d + l, nx)] += -2*beta*x[(N+1)*d + i*d + k]*x[(N+1)*d + i*d + l];
			}	
		}
	}
}




void GM(double *res, double *x, int N, int d, int nx)
{
	int i, j, k, l, s;

	for(i=0; i<N+1; i++)
	{
		for(j=0; j<N+1; j++)
		{
			if(i!=j)
			{
				for(k=0; k<d; k++)
				{
					for(l=0; l<d; l++)
					{
						res[map((N+1)*d + i*d + k, j*d + l, nx)] += -dmdy(x, i, j, k, l, N, d, Ca, la, Cr, lr)/(N+1);
					}
				}
			}
			else
			{
				for(k=0; k<d; k++)
				{
					for(l=0; l<d; l++)
					{
						double temp=0;
						for(s=0; s<N+1; s++)
						{
							if(s!=i)
							{
								temp += dmdx(x, i, s, k, l, N, d, Ca, la, Cr, lr);
							}
						}
						res[map((N+1)*d + i*d + k, j*d + l, nx)] += -temp/(N+1);
					}
				}

			}
		}
	}
}

void GL(double *res, double *x, int N, int d, int nx)
{
	int i, k, l;

	for(i=1; i<N+1; i++)
	{
		for(k=0; k<d; k++)
		{
			for(l=0; l<d; l++)
			{
				res[map((N+1)*d + i*d + k, i*d + l, nx)] += -gamma1*dmdx(x, i, 0, k, l, N, d, Ca0, la0, Cr0, lr0);
				res[map((N+1)*d + i*d + k, l, nx)] += -gamma1*dmdy(x, i, 0, k, l, N, d, Ca0, la0, Cr0, lr0);

			}
		}
	}
}


double dmdx(double *x, int i, int j, int k, int l, int N, int d, double Ca, double la, double Cr, double lr)// i,j = 1...N+1;    k,l = 1...d;
{
	int z;

	double temp[d];
	for(z=0; z<d; z++)
	{	
		temp[z] = x[i*d+z] - x[j*d+z];
	}
	double ro = norm(temp, d);

	

	double res = -(x[i*d+l] - x[j*d+l])*(x[i*d+k] - x[j*d+k])/(ro*ro*ro);
	if(k==l)
	{
		res += 1/ro;
	}

	double m1 = (Ca/la*exp(-ro/la) - Cr/lr*exp(-ro/lr));
	res *= m1;

	res += (x[i*d+l] - x[j*d+l])*(x[i*d+k] - x[j*d+k])*(Cr/(lr*lr)*exp(-ro/lr) - Ca/(la*la)*exp(-ro/la))/(ro*ro);

	return res;
}
double dmdy(double *x, int i, int j, int k, int l, int N, int d, double Ca, double la, double Cr, double lr)// i,j = 1...N+1;    k,l = 1...d;
{
	return - dmdx(x, i, j, k, l, N, d, Ca, la, Cr, lr);
}





// Objective /////////////////////////////////////////////////////////////////////////////////////////////////////////////

double l1(double *x, int N, int d)
{
	int i, k;

	double temp = 0;
	for(i=1; i<N+1; i++)
	{
		double temp1[2];
		for(k=0; k<d; k++)
		{
			temp1[k] = x[k] - x[i*d + k];
		}
		double temp3 = norm(temp1, d);

		temp += temp3*temp3*temp3*temp3;
	}

	return 0.5*mu*temp;
}

double Gl1(double *res,  double *x, int N, int d, int nx)
{
	 int j, k;

	 double temp[2];
	 for(j=1; j<N+1; j++)
	 { 

		 double temp1[2];
		 for(k=0; k<d; k++)
		 {
	 		temp1[k] = x[k] - x[j*d + k];
         }
		 double temp3 = norm(temp1, d);

		 for(k=0; k<d; k++)
		 {
	         res[map(nx-1, j*d + k, nx)] = -2*mu*temp3*temp3*temp1[k];
	         

	         temp[k] += temp3*temp3*temp1[k];
		 }
	 }

	 for(k=0; k<d; k++)
	 {
	 	res[map(nx-1, k, nx)] = 2*mu*temp[k];
	 }	
}




//  Destination /////////////////////////////////////////////////////////////////////////////////////////////////////////

// free memory after use!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double xdes(double t, int k)
{
	// double *res = (double*) malloc(d*sizeof(double));
	// // double res[d];

	// int k;
	// for(k=0; k<d; k++)
	// {
	// 	res[k] = t;
	// }
	// return res;
	if(k == 0) return  -t;
	else if (k == 1) return t;
	
	return t;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////

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
