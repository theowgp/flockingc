#include "plugin1.h"

void pad(double *f, int n)
{
     // printf("pad = \n");
     int i;
     for(i=0; i<n; i++)
     {
       printf("%f  ", f[i]);
       // printf("%f  ", *f++);
     } 
     printf("\n");
}

void pmd(double *f, int n, int m)
{
     int i, j;
     // printf("pmd = \n");
     for(i=0; i<n; i++)
     {
     	printf("\n");
     	for(j=0; j<m; j++)
     	{
     		printf("%f    ", f[map(i, j, n)]);
       
     	}
       
     } 
     printf("\n");
}

void setmt0(double *f, int n, int m)
{
     int i, j;
     for(i=0; i<n; i++)
	 {
		 for (j=0; j<m; j++)
		 {
			 f[map(i, j, n)] = 0;
		 }
	 }
}