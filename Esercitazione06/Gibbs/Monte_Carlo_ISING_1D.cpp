/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization

  for(int k = 0; k < 2; k++){

    //ciclo per variare il valore del campo magnetico (h)
    h = 0.0 + k*0.02;
    cout << "External field = " << h << endl << endl;

    //setto i parametri del ciclo per valutare l'andamento delle grandezze al variare della temperatura
    int npoint = 20;
    double tmin = 0.5;
    double tmax = 2.0;

    for(int i = 0; i <= npoint; i++){

      temp = tmin + i*(tmax - tmin)/(double)npoint;
      beta = 1.0/temp;
      cout << "Temperature = " << temp << endl << endl;
      m_real = MagReal();

      //comincio la simulazione con la nuova temperatura settata
      nstart = 1; //tengo il conto del numero di (re)start
      Config(); //configuro il sistema nello stato iniziale
      Equilibrium(); //fase di equilibrazione
      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move(metro);
          Measure();
          Accumulate(); //Update block averages
        }
        Averages(iblk);   //Print results for current block
      }
      ConfFinal(); //Write final configuration
    }
  }
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  //cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> read;  // if=1 Read from file "config.0"

  ReadInput >> nstep_eq;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables
}

void Config(void)
{
  if(read == 0){

    //initial configuration
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 1) s[i] = 1;
      else s[i] = -1;
    }

  }else{

    ifstream ReadConfig("config.0");

    //reading initial configuration
    for (int i=0; i<nspin; ++i)  ReadConfig >> s[i];
  }
}

void Equilibrium(){
  double m, err;
  do{
      double sumM = 0.0, sumM2 = 0.0;
      for(int istep = 1; istep <= nstep_eq; ++istep){
         Move(metro);
         m = MagIst();
         sumM += m;
         sumM2 += m*m;
      }

      sumM /= nstep_eq;
      sumM2 /= nstep_eq;

      m = sumM;
      err = sqrt((sumM2 - sumM*sumM)/nstep_eq);

      //confronto le fluttuazioni della temperatura istantanea con quella richiesta
      if(abs(m - m_real) < err) restart = false;
      else restart = true;

      ConfFinal();
      read = true;
      OverWrite(); //Sovrascrivo le coordinate iniziali con quelle finali
      nstart ++;
      Config(); //Inizialization
    }while(restart == true);
  cout << "Equilibrium state reached..." << endl;
  cout << "Number of simulations to reach the equilibrium state: " << (nstart-1)*nstep_eq << endl << endl;
}

void Move(int metro)
{
  int o;
  double energy_old, energy_new, sm, deltaE;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      sm = s[o];
      energy_old = Boltzmann(sm, o);
      energy_new = Boltzmann(-1*sm, o);
      deltaE = energy_new - energy_old;
      double p = exp(-beta*deltaE);

      if(deltaE <= 0){

        s[o] = sm*(-1);
        accepted++;

      }else{

        double r = rnd.Rannyu();

        if(r <= p){

          s[o] = sm*(-1);
          accepted++;
        }
      }

    }
    else //Gibbs sampling
    {
      sm = s[o];
      energy_old = Boltzmann(sm, o);
      energy_new = Boltzmann(-1*sm, o);
      deltaE = energy_new - energy_old;
      double p_up = 1./(1. + exp(-sm*beta*deltaE));
      double r = rnd.Rannyu();

      if(r <= p_up){

        s[o] = 1;
        accepted++;

      }else{

        s[o] = -1;
        accepted++;

      }
    }

    attempted++;

  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    blk_av[iu] /= blk_norm;
    stima_u = blk_av[iu]/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u = Error(glob_av[iu],glob_av2[iu],iblk);

    blk_av[ic] /= blk_norm;
    stima_c = beta*beta*(blk_av[ic] - blk_av[iu]*blk_av[iu]);
    stima_c /= (double)nspin; //Heat Capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c = Error(glob_av[ic],glob_av2[ic],iblk);

    blk_av[im] /= blk_norm;
    stima_m = blk_av[im]/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m = Error(glob_av[im],glob_av2[im],iblk);

    blk_av[ix] /= blk_norm;
    stima_x = beta*(blk_av[ix] - blk_av[im]*blk_av[im]);
    stima_x /= (double)nspin; //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x = Error(glob_av[ix],glob_av2[ix],iblk);

    //carico su file
    if(h == 0){

      Ene.open("output.ene.0",ios::app);
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();

      Heat.open("output.heat.0",ios::app);
      Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Heat.close();

      Chi.open("output.chi.0",ios::app);
      Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();

    }else{

      Mag.open("output.mag.0",ios::app);
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();

    }

    if(iblk == nblk){

      ofstream Ufinal, Cfinal, Mfinal, Xfinal;
      ofstream Acc;

      double u,c,m,x,a;

      u = glob_av[iu]/(double)iblk;
      c = glob_av[ic]/(double)iblk;
      m = glob_av[im]/(double)iblk;
      x = glob_av[ix]/(double)iblk;
      a = (double)accepted/attempted;

      //carico su file
      if(h == 0){

        Ufinal.open("ave_ene.final",ios::app);
        Ufinal << temp << " " << u << " " << err_u << endl;
        Ene.close();

        Cfinal.open("ave_heat.final",ios::app);
        Cfinal << temp << " " << c << " " << err_c << endl;
        Cfinal.close();

        Xfinal.open("ave_chi.final",ios::app);
        Xfinal << temp << " " << x << " " << err_x << endl;
        Xfinal.close();

        Acc.open("acc.final",ios::app);
        Acc << temp << " " << a << endl;
        Acc.close();

      }else{

        Mfinal.open("ave_mag.final",ios::app);
        Mfinal << temp << " " << m << " " << err_m << endl;
        Mfinal.close();

      }
  
      cout << "Acceptance rate = " << a << endl;
      cout << "Inernal energy  = " << u << endl;
      cout << "Heat Capacity   = " << c << endl;
      cout << "Magnetization   = " << m << endl;
      cout << "Susceptibility  = " << x << endl;

      cout << "----------------------------" << endl << endl;

    }

}


void ConfFinal(void)
{
  ofstream WriteConf;

  //cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk == 1)

      return 0;

    else

      return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void OverWrite(){
  ofstream WriteConf;
  ifstream ReadConf;

  double s;

  WriteConf.open("config.0");
  ReadConf.open("config.final");

  for (int i=0; i<nspin; ++i){
    ReadConf >> s;
    WriteConf << s << endl;
  }

  WriteConf.close();
  ReadConf.close();

  return;
}

double MagIst(){

  double m = 0.0;

  for (int i=0; i<nspin; ++i) m += s[i];

  return m/nspin;
}

double MagReal(){

  double beta = 1./temp;
  double J = 1.;
  double nspin = 50;
  double h = 0.02;

  double l1 = exp(beta*J)*cosh(beta*h) + sqrt(exp(2.*beta*J)*cosh(beta*h)*cosh(beta*h) - 2.*sinh(2.*beta*J));
  double l2 = exp(beta*J)*cosh(beta*h) - sqrt(exp(2.*beta*J)*cosh(beta*h)*cosh(beta*h) - 2.*sinh(2.*beta*J));
  double Z = pow(l1,nspin) + pow(l2,nspin);
  double k1 = pow(l1,nspin-1)*(1. + (exp(beta*J)*cosh(beta*h))/sqrt(exp(2.*beta*J)*cosh(beta*h)*cosh(beta*h) - 2.*sinh(2.*beta*J)));
  double k2 = pow(l2,nspin-1)*(1. - (exp(beta*J)*cosh(beta*h))/sqrt(exp(2.*beta*J)*cosh(beta*h)*cosh(beta*h) - 2.*sinh(2.*beta*J)));

  double M = (exp(beta*J)*sinh(beta*h)/(Z))*(k1 + k2);

  return M;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
