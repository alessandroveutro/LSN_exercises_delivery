#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "random.h"
#include "funzioni.h"

using namespace std;

double Error(double *av, double *av2, int N);

void UniformSampling(double *x, Random rnd, int dim);
void ImportanceSampling(double *x, Random rnd, int dim);

void BlockMethod(double *mean, double *err_mean, double *x, Function1D *f, int n_blocks, int L);

