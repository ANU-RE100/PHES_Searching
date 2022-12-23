#ifndef FSOLVE
#define FSOLVE
#include <bits/stdc++.h>

void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
  double delta, double x[], double wa1[], double wa2[] );
double enorm ( int n, double x[] );
void fdjac1 ( std::function<void (int, double *, double *)> fcn, 
  int n, double x[], double fvec[], double fjac[], int ldfjac,
  int ml, int mu, double epsfcn, double wa1[], double wa2[] );
int fsolve ( std::function<void (int, double *, double *)> fcn, int n, 
  double x[], double fvec[], double tol, double wa[], int lwa );
int hybrd ( std::function<void (int, double *, double *)> fcn, 
  int n, double x[], double fvec[], double xtol, int maxfev, int ml, 
  int mu, double epsfcn, double diag[], int mode, double factor, 
  int nfev, double fjac[], int ldfjac, double r[], int lr, double qtf[], 
  double wa1[], double wa2[], double wa3[], double wa4[] );
void qform ( int m, int n, double q[], int ldq );
void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
  int lipvt, double rdiag[], double acnorm[] );
void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] );
bool r1updt ( int m, int n, double s[], int ls, double u[], double v[], double w[] );

#endif
