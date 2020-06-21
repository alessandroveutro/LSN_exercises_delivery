#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include "random.h"
#include "point3D.h"

using namespace std;

double Error(double *av, double *av2, int N);

double Psi_100(Point3D x);
double Psi_210(Point3D x);

void MoveUnif_100(Point3D& point, Random& rnd, int& accepted);
void MoveUnif_210(Point3D& point, Random& rnd, int& accepted);
void MoveGauss_100(Point3D& point, Random& rnd, int& accepted);
void MoveGauss_210(Point3D& point, Random& rnd, int& accepted);

void EquilUnif_100(int nstep, Point3D point, double r_real, Random rnd);
void EquilUnif_210(int nstep, Point3D point, double r_real, Random rnd);
void EquilGauss_100(int nstep, Point3D point, double r_real, Random rnd);
void EquilGauss_210(int nstep, Point3D point, double r_real, Random rnd);

void DataUnif_100(int n_throws, Point3D point, double *r, Random rnd);
void DataUnif_210(int n_throws, Point3D point, double *r, Random rnd);
void DataGauss_100(int n_throws, Point3D point, double *r, Random rnd);
void DataGauss_210(int n_throws, Point3D point, double *r, Random rnd);

void BlockMethod(double *mean, double *err_mean, double *data, int n_blocks, int L);
