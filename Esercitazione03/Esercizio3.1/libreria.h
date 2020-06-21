#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "random.h"

using namespace std;

double Error(double *av, double *av2, int N);

double N(double x);
double C_BlackScholes(double S0, double K, double T, double r, double sigma);
double P_BlackScholes(double S0, double K, double T, double r, double sigma);

void FillRandomGauss(double *z, Random rnd, int dim, double mean, double sigma);
void BlockMethod(double *mean, double *err_mean, double *x, int n_blocks, int L);

