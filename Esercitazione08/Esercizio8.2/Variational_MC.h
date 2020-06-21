#ifndef __variabili__
#define __variabili__

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props = 1000;
int ih,ip,n_props;
double mu,sigma,x;
double walker[m_props];

// simulation
int nstep, nblk;
double delta;
bool opt;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_ene,err_ene,stima_prob,err_prob;

//histogram
int nbins = 80;
double xmin = -5.;
double xmax = 5.;
double step = (xmax-xmin)/(double)nbins;
//double histo[nbins];

//functions
void Input(void);
void Equilibrium(void);
void Variational(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
void WriteAve(int);
double Hamilton(double);
double Psi_Trial(double);
double Psi_Trial_Second(double);
double V(double);
double Error(double,double,int);

#endif
