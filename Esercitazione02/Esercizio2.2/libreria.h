#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include "random.h"
#include "randwalk.h"

using namespace std;

double Error(double *av, double *av2, int N);

void FillRandomDiscrete(double *r, Random random, int n_steps, int n_walks);
void FillRandomContinuum(double *r1, double *r2, Random random, int n_steps, int n_walks);

void RW_Discrete(double *mean, double *err_mean, RandomWalk *rw, double *rand, double a, int n_steps, int n_walks);
void RW_Continuum(double *mean, double *err_mean, RandomWalk *rw, double *rand1, double *rand2, double a, int n_steps, int n_walks);
