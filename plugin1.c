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
     		printf("%f    ", f[map(i, j, m)]);
       
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
			 f[map(i, j, m)] = 0;
		 }
	 }
}

void wtf(double *state, int n, int ns, int nx, int N, int d)
{


      FILE *f = fopen("state.txt", "w");
      if (f == NULL)
      {
          printf("Error opening file!\n");
        
      }
 
      
      
      fprintf(f, " %d\n", n);
      fprintf(f, " %d\n", ns);
      fprintf(f, " %d\n", nx);
      fprintf(f, " %d\n", N);
      fprintf(f, " %d\n", d);
 
      /* print integers and floats */
      int i;
      for(i=0; i<n*ns*nx+nx; i++)
      {
          fprintf(f, "%f\n", state[i]);    
      }
      
      fclose(f);
}