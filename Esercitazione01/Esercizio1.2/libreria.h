#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "random.h"

using namespace std;

void StandardDice(double *r, Random rnd, int dim);
void ExpDice(double *r, Random rnd, int dim, double lambda);
void CauchyDice(double *r, Random rnd, int dim, double mu, double gamma);
