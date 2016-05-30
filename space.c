#include "space.h"

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

double normL2(double *v, int n, double h)
{
  double temp = 0;

  int i;
  for(i = 0; i<n; i++)
  {
    temp += h*v[i]*v[i];
  } 

  return sqrt(temp);
}