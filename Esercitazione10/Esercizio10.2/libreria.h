#include "path.h"

int Select(int N, Random *rnd);
int PBC(int k, int N);
void Crossover(Path *path, int i, int N, int l, int n);
void Cross(Path *path, int *genes, int i1, int i2, int l, int n);
void Evolution(int N, int l, Path *path, Random *rnd, map<int,City> mapping);
void Order(int N, Path *path);
void Equal(Path *path1, Path *path2, int i1, int i2, map<int,City> mapping);
void WriteBest(Path path, int istep, int rank);
void WriteMean(Path *path, int istep, int rank);
void WriteFinal(Path path, int l, map<int,City> mapping, int rank);
