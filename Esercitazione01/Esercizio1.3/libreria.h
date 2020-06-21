#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "random.h"
#include "needle.h"

using namespace std;

double Error(double *av, double *av2, int N);
void ThrowNeedle(Needle *needle, Random rnd, double d, double l, int n_throws);
void BlockMethod(double *mean, double *err_mean, Needle *needle, int n_blocks, int L, double d, double l);
