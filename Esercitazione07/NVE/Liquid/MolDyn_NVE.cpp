/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <ostream>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  nstart = 1; //tengo il conto del numero di (re)start
  Input(); //carico i dati da input.dat
  Read(); //leggo le configurazioni delle particelle da file
  Equilibrium(); //Fase di equilibrazione

  int nconf = 1;

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(); //Move particles with Verlet algorithm
      Measure(istep); //Properties measurement
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  if(old == true) OldFinal(); //Write pre-final configuration to restart

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Input
void Input(void){ //Prepare all stuff for the simulation

   ifstream ReadInput;
   ReadInput.open("input.dat"); //Read input

   cout << "Classic Lennard-Jones fluid        " << endl;
   cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
   cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
   cout << "The program uses Lennard-Jones units " << endl;

   ReadInput >> temp;
   cout << "Desired temperature = " << temp << endl;
   ReadInput >> npart;
   cout << "Number of particles = " << npart << endl;
   ReadInput >> rho;
   cout << "Density of particles = " << rho << endl;
   vol = (double)npart/rho;
   cout << "Volume of the simulation box = " << vol << endl;
   box = pow(vol,1.0/3.0);
   cout << "Edge of the simulation box = " << box << endl;

   ReadInput >> rcut;
   ReadInput >> delta;
   ReadInput >> nstep;
   ReadInput >> iprint;

   cout << "The program integrates Newton equations with the Verlet method " << endl;
   cout << "Time step = " << delta << endl;
   cout << "Number of steps (in each block) = " << nstep << endl;

   ReadInput >> nstep_eq;
   int n;
   ReadInput >> n;
   if(n == 1) old = true;
   else old = false;

   ReadInput >> nblk;

   cout << "Number of blocks = " << nblk << endl << endl;

   ReadInput.close();

   seed = 1;    //Set seed for random numbers
   srand(seed); //Initialize random number generator

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //Virial (Pressure)
  n_props = 5; //Number of observables

//measurement of g(r)
  igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

   return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Read
void Read(){
  ifstream ReadConf;

  //Read initial configuration
  if(nstart == 1)
    cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
//*************************************************************************************************************************************************
  if(old == true){
    ifstream ReadOld;
    if(nstart == 1)
      cout << "Read initial configuration from file old.0 " << endl << endl;
    ReadOld.open("old.0");
    for (int i=0; i<npart; ++i){
      ReadOld >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadOld.close();
    Move(); //eseguo un passo con l'algoritmo di Verlet per ottenere r(t+dt)
    double sumv2 = 0.0, fs;
    double v2;

    for (int i=0; i<npart; ++i){
      vx[i] = Pbc(x[i] - xold[i])/delta;
      vy[i] = Pbc(y[i] - yold[i])/delta;
      vz[i] = Pbc(z[i] - zold[i])/delta;

      v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      sumv2 += v2;
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
//*************************************************************************************************************************************************
  }else{
   //Prepare initial velocities
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

    Move(); //eseguo un passo con l'algoritmo di Verlet per ottenere r(t+dt)
    sumv2 = 0.0;
    double v2;

    for (int i=0; i<npart; ++i){
      vx[i] = Pbc(x[i] - xold[i])/delta;
      vy[i] = Pbc(y[i] - yold[i])/delta;
      vz[i] = Pbc(z[i] - zold[i])/delta;

      v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      sumv2 += v2;
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }

   return;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Equilibrium
void Equilibrium(){
  double t,u;
  for(int i = 0; i < 6; i++){
      ofstream WriteT, WriteU;
      WriteT.open("equil_temp.out",ios::app);
      WriteU.open("equil_epot.out",ios::app);
      for(int istep = 1; istep <= nstep_eq; ++istep){
         //Scrivo su file la temperatura istantanea
         t = Temperature();
         u = Potential();

         WriteT << t << endl;
         WriteU << u << endl;

         Move();

      }

      WriteT.close();
      WriteU.close();
      ConfFinal();
      old = true;
      OldFinal();
      OverWrite(); //Sovrascrivo le coordinate iniziali con quelle finali
      nstart ++;
      Read(); //Inizialization
    }
  cout << "Equilibrium state reached..." << endl;
  cout << "Number of simulations to reach the equilibrium state: " << (nstart-1)*nstep_eq << endl << endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Move
void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Temperature
double Temperature(){
  double t = 0.0;
  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  return (2.0 / 3.0) * t/(double)npart;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Potential
double Potential(){
  double v = 0.0, vij;
  double dx, dy, dz, dr;
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

  v = 4.0 * v; //Potential energy per particle
  v /= (double)npart;

  return v;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Force
double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Measure
void Measure(int istep){ //Properties measurement
  double v, t, w, vij, wij;
  double dx, dy, dz, dr;
  double rmin, rmax, dVol;

  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);

  //reset observables
  v = 0.0; 
  t = 0.0;
  w = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
     for(int k = 0; k < nbins; k++){

       rmin = k*bin_size;
       rmax = rmin + bin_size;

       if(dr >= rmin and dr < rmax)

         walker[igofr+k] += 2.0;

     }

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

//Potential energy
       v += vij;
//Virial
       w += wij;
     }
    }          
  }

  v = 4.0 * v; //Potential energy per particle
  w = 48.0 * w / 3.0; //Virial

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = v/(double)npart; //Potential energy per particle
  walker[ik] = t/(double)npart; //Kinetic energy per particle
  walker[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
  walker[ie] = (t+v)/(double)npart; //Total energy per particle
  walker[iw] = w; //Virial

  for(int i = 0; i < nbins; i++){

    rmin = i*bin_size;
    rmax = rmin + bin_size;   
    dVol = (4.*pi/3.)*(pow(rmax,3) - pow(rmin,3));
    walker[igofr+i] /= (rho*npart*dVol);

  }

  if(istep%100 == 0){

    double p = rho*temp + (1.0/vol)*w;

    Epot << walker[iv]  << endl;
    Ekin << walker[ik]  << endl;
    Temp << walker[it] << endl;
    Etot << setprecision(8) << walker[ie] << endl;
    Pres << p << endl;

  }

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Reset
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
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Accumulate
void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Averages
void Averages(int iblk) //Print results for current block
{
   double r;
   ofstream Epot, Ekin, Etot, Temp, Pres;
   ofstream Gofr, Gave;
    
    cout << "Block number " << iblk << endl;

    Epot.open("ave_epot.out",ios::app);
    Ekin.open("ave_ekin.out",ios::app);
    Temp.open("ave_temp.out",ios::app);
    Etot.open("ave_etot.out",ios::app);
    Pres.open("ave_pres.out",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);

    stima_pot = blk_av[iv]/blk_norm; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_epot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_kin = blk_av[ik]/blk_norm; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_ekin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    stima_pres = rho * temp + (blk_av[iw]/blk_norm)/ vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_pres=Error(glob_av[iw],glob_av2[iw],iblk);

    for(int i = igofr; i < igofr+nbins; i++){ //g(r) histogramm

      stima_gdir = blk_av[i]/(double)blk_norm;
      glob_av[i] += stima_gdir;
      glob_av2[i] += stima_gdir*stima_gdir;

    }

    Epot << iblk <<  " " << stima_pot << " " << glob_av[iv]/(double)iblk << " " << err_epot << endl;
    Ekin << iblk <<  " " << stima_kin << " " << glob_av[ik]/(double)iblk << " " << err_ekin << endl;
    Etot << iblk <<  " " << stima_etot << " " << setprecision(8) << glob_av[ie]/(double)iblk << " " << err_etot << endl;
    Temp << iblk <<  " " << stima_temp << " " << glob_av[it]/(double)iblk << " " << err_temp << endl;
    Pres << iblk <<  " " << stima_pres << " " << glob_av[iw]/(double)iblk << " " << err_pres << endl;

    //g(r)
    //mean in each block
    Gofr << iblk;
    for(int i = igofr; i < igofr+nbins; i++){

      Gofr << " " << stima_gdir;

    }

    Gofr << endl;

    //final average
    if(iblk == nblk){

      int k = 0;

      for(int i = igofr; i < igofr+nbins; i++){

        r = k*bin_size;
        err_gdir = Error(glob_av[i], glob_av2[i], nblk);
        Gave << r << " " << glob_av[i]/(double)nblk << " " << err_gdir << endl;
        k++;

      }      
    }

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
    Gofr.close();
    Gave.close();

    cout << "----------------------------" << endl << endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  ConfFinal
void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << endl << "Print final configuration to file config.final " << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  OldFinal
void OldFinal(void){ //Write final configuration
  ofstream WriteOld;

  cout << "Print pre-final configuration to file old.final " << endl << endl;
  WriteOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteOld.close();
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  OverWrite
void OverWrite(){
  ofstream WriteConf, WriteOld;
  ifstream ReadConf, ReadOld;

  WriteConf.open("config.0");
  ReadConf.open("config.final");

  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    WriteConf << x[i] << "   " <<  y[i] << "   " << z[i] << endl;
  }

  WriteConf.close();
  ReadConf.close();

  if(old == true){
    ofstream WriteOld;
    ifstream ReadOld;

    WriteOld.open("old.0");
    ReadOld.open("old.final");

    for (int i=0; i<npart; ++i){
      ReadOld >> xold[i] >> yold[i] >> zold[i];
      WriteOld << xold[i] << "   " <<  yold[i] << "   " << zold[i] << endl;
    }

    WriteOld.close();
    ReadOld.close();
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  ConfXYZ
void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Pbc
double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
  return r - box * rint(r/box);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  Error
double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt(abs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
