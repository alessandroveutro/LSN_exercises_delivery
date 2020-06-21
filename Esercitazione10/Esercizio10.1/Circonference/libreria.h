#include "path.h"

int PBC(int k, int N);
double Boltzmann(double temp, double Lold, double Lnew);

void Annealing(int l, Path *path, Random *rnd, map<int,City> mapping, double temp);

void Equal(Path *path1, Path *path2, map<int,City> mapping);
void WriteBest(Path *path, double temp);
void WriteFinal(Path *path, int l, map<int,City> mapping);
