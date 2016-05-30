#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#include "space.h"


//  f ////////////////////////////////////////////////////////////////////////////////////////////////////////
void S(double *res, double *v, int N, int d);
void Si(double *res, double *v, int d);

void fx(double *res, double *v, int N, int d);




//  Gf ///////////////////////////////////////////////////////////////////////////////////////////////////////
void GS(double *res, double *x, int N, int d, int nx);

void GM(double *res, double *x, int N, int d, int nx);
double dmdx(double *x, int i, int j, int k, int l, int N, int d, double Ca, double la, double Cr, double lr);
double dmdy(double *x, int i, int j, int k, int l, int N, int d, double Ca, double la, double Cr, double lr);


// Objective /////////////////////////////////////////////////////////////////////////////////////////
double l1(double *x, int N, int d);

 //  Destination /////////////////////////////////////////////////////////////////////////////////////////////////////////
double xdes(double t, int k);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int map(int i, int j, int nx);
int* imap(int k, int nx);

