/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props = 105;
int n_props;
int iv,ik,it,ie,iw,igofr;
double bin_size,nbins;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres,stima_gdir;
double walker[m_props];

//configuration
const int m_part = 108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nstep_eq, nstart;
double delta;
bool old;
//bool restart;

//pigreco
const double pi=3.1415927;

// averages
int nblk;
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double err_epot, err_ekin, err_etot, err_temp, err_pres, err_gdir;

//functions
void Input(void);
void Read(void);
void Equilibrium(void);
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void OldFinal(void);
void OverWrite(void);
void ConfXYZ(int);
void Measure(int);
double Temperature(void);
double Potential(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
