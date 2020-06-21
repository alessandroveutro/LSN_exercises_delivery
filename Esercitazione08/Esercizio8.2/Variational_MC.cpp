#include "Variational_MC.h"


int main(int argc, char *argv[]){

	Input(); //Inizialization
	x = 0.0; //punto iniziale
	Variational();
	Equilibrium();

	for(int iblk = 1; iblk <= nblk; ++iblk){

			Reset(iblk);   //Reset block averages

			for(int istep = 1; istep <= nstep; ++istep){

				Move();
				Measure();
				Accumulate(); //Update block averages

			}

			Averages(iblk);   //Print results for current block
			WriteAve(iblk);

	}

	//plot grafico della distribuzione di probabilità di prova minimizzata
	for (int i = 0; i < 100; i++){

		double k = -3 + (double)i*6./100.;
		double psi = Psi_Trial(k);

	}

	rnd.SaveSeed();

	return 0;

}

void Input(void){

	ifstream ReadInput;

	cout << "Variational Method                 " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Potential 1D v(x) = x^4 - (5/2)*x^2" << endl << endl;
	cout << "The program uses classical units " << endl;

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

	ReadInput >> mu;
	ReadInput >> sigma;
	ReadInput >> delta;
	ReadInput >> nblk;
	ReadInput >> nstep;

	cout << "The program perform Metropolis moves with uniform translations" << endl;
	cout << "Moves parameter = " << delta << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;

	ReadInput.close();

	n_props = 101;
	ih = 0;
	ip = 1;

}

void Equilibrium(void){

	if(opt) cout << "Equilibration Phase starting..." << endl;

	for (int i = 0; i < 1000; i++) Move();

	if(opt){

		cout << "Equilibrium reached!" << endl;
		cout << "----------------------------" << endl << endl;

	}
}

void Variational(void){

	double H_blk;
	double H_test = 0.0,mu_test = mu,sigma_test = sigma;
	int N = 100;

	opt = false;

	double mu_min = 0.7;
	double mu_max = 0.8;
	double sigma_min = 0.55;
	double sigma_max = 0.65;

	for(int i = 0; i < N; i++){

		mu = mu_min + i*(mu_max-mu_min)/N;
		sigma = sigma_min + i*(sigma_max-sigma_min)/N;
		Equilibrium();

		//Simulation
		for(int iblk = 1; iblk <= nblk; ++iblk){

			Reset(iblk);   //Reset block averages

			for(int istep = 1; istep <= nstep; ++istep){

				Move();
				Measure();
				Accumulate(); //Update block averages

			}

			Averages(iblk);   //Print results for current block

		}

		H_blk = glob_av[ih]/(double)nblk;

		if(i == 0) H_test = H_blk;
		else{

			if(H_blk < H_test){

				H_test = H_blk;
				mu_test = mu;
				sigma_test = sigma;

			}
		}
	}

	opt = true;
	//carico il valore dei parametri ottimizzati
	mu = mu_test;
	sigma = sigma_test;

	cout << "Trial function minimized!!" << endl;
	cout << "Optimized parameters found!!" << endl;
	cout << "\t mu    = " << mu << endl;
	cout << "\t sigma = " << sigma << endl;
	cout << "----------------------------" << endl << endl;

}

void Reset(int iblk){
	
	if(iblk == 1){

		for(int i = 0; i < n_props; ++i){

			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i = 0; i < n_props; ++i){

		blk_av[i] = 0;
	}

	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}

void Accumulate(void){

	for(int i=0; i < n_props; ++i){

		blk_av[i] = blk_av[i] + walker[i];

	}

	blk_norm = blk_norm + 1.0;
}

void Averages(int iblk){
    
	if(opt){

		cout << "Block number " << iblk << "/" << nblk << endl;
		cout << "Acceptance rate " << accepted/attempted << endl << endl;
		cout << "----------------------------" << endl << endl;

	}

	stima_ene = blk_av[ih]/blk_norm;
	glob_av[ih] += stima_ene;
	glob_av2[ih] += stima_ene*stima_ene;
	err_ene = Error(glob_av[ih],glob_av2[ih],iblk);

	for(int i = ip; i < ip+nbins; i++){ //prob(r) histogramm

		stima_prob = blk_av[i]/(double)blk_norm;
		glob_av[i] += stima_prob;
		glob_av2[i] += stima_prob*stima_prob;

	}
}

void Move(void){

	double xold = x;
	double xnew = rnd.Rannyu(x-delta,x+delta);

	double p_new = pow(Psi_Trial(xnew),2.0);
	double p_old = pow(Psi_Trial(xold),2.0);

	double t = p_new/p_old;
	double alpha;
	
	if(t > 1)
		
		alpha = 1;

	else

		alpha = t;

	double s = rnd.Rannyu();

	if(s < alpha){

		x = xnew;
		accepted++;

	}

	attempted++;
}

void Measure(void){

	walker[ih] = Hamilton(x);

	for(int k = ip; k < nbins+ip; k++) walker[k]=0.0;

	for(int k = 0; k < nbins; k++)

		if(x >= xmin + k*step and x < xmin + (k+1)*step)

			walker[ip+k] += 1;

}

void WriteAve(int iblk){

	//Print results for current block
	ofstream Eave,Pave;
	Eave.open("average.ene.0",ios::app);
	Pave.open("average.prob.0",ios::app);

	//Scrivo su file il valor medio dell'Hamiltoniana ad ogni blocco
	Eave << iblk <<  " " << stima_ene << " " << glob_av[ih]/(double)iblk << " " << err_ene << endl;

	//final average
	if(iblk == nblk){

		int k = 0;

		for(int i = ip; i < ip+nbins; i++){

			x = xmin + k*step;
			err_prob = Error(glob_av[i], glob_av2[i], nblk);
			Pave << x << " " << glob_av[i]/(double)nblk/step << " " << err_prob/step << endl;
			k++;

		}      
	}

	Eave.close();
	Pave.close();

}

double Hamilton(double y){

	double psi = Psi_Trial(y);
	double pot = V(y)*psi;
	double kin = -0.5*Psi_Trial_Second(y);
	double H = (kin + pot)/psi;

	return H;
}

double Psi_Trial(double y){

	//parametri liberi mu e sigma
	double a1 = -(y - mu)*(y - mu)/(2.0*sigma*sigma);
	double a2 = -(y + mu)*(y + mu)/(2.0*sigma*sigma);
	double psi = exp(a1) + exp(a2);

	return psi;

}

double Psi_Trial_Second(double y){

	//parametri liberi mu e sigma
	double a1 = -(y - mu)*(y - mu)/(2.0*sigma*sigma);
	double a2 = -(y + mu)*(y + mu)/(2.0*sigma*sigma);
	double b1 = (y - mu)*(y - mu)/pow(sigma,4.0);
	double b2 = (y + mu)*(y + mu)/pow(sigma,4.0);

	return b1*exp(a1) + b2*exp(a2) - (1.0/(sigma*sigma))*Psi_Trial(y);

}

double V(double y){

	return pow(y,4.0) - 2.5*pow(y,2.0);

}

double Error(double sum, double sum2, int iblk){

    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));

}
