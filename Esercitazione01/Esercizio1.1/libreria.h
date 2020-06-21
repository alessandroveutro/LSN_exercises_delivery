#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "random.h"

using namespace std;

double Error(double *av, double *av2, int N);
void FillRandom(double *r, Random rnd, int dim, double min, double max);
void BlockMean(double *mean, double *err_mean, double *r, int N, int L);
void BlockVariance(double *var, double *err_var, double *r, int N, int L);
double ChiQuadro(int *n, double *r, int n_step, int n_throws, int l);


