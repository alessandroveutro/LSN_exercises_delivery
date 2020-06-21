#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

double Mean(double *x, int n);
double Variance(double *x, int n);
double Autocorr(double *x, int n, int t);

double Error(double ave, double ave2, int n);
double ErrorBM(double *x, int nblk, int L);


