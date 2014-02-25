/* File foo.c */
#include <complex.h>
#include <math.h>

double _Complex foo(int n, double a, double *V, double _Complex *expon, double _Complex *calcdata) {
	int i;
	for (i=0;i<n;i++) {
		calcdata[i] = V[i]*cexp(a*expon[i]);
  	}
}