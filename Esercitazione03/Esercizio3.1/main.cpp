#include "random.h"
#include "libreria.h"

int main(int argc, char *argv[]){

	//setto il generatore di numeri casuali
	Random rnd; 
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	if (Primes.is_open()){
		
		Primes >> p1 >> p2 ;
	
	} else cerr << "PROBLEM: Unable to open Primes" << endl;

	Primes.close();

	ifstream input("seed.in");
	string property;

	if (input.is_open()){
	while ( !input.eof() ){
		input >> property;
		if( property == "RANDOMSEED" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
		}
	}

	input.close();

	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	//inizio esercizio

	int n_throws = 10000;
	int n_blocks = 100;
	int L = n_throws/n_blocks;

	double S0 = 100;
	double T = 1;
	double K = 100;
	double r = 0.1;
	double sigma = 0.25;

	//Calcolo i valori attesi di call e put, utilizzando il modello Black-Scholes

	double C_BS = C_BlackScholes(S0, K, T, r, sigma);
	double P_BS = P_BlackScholes(S0, K, T, r, sigma);

	double *C0 = new double[n_throws];
	double *P0 = new double[n_throws];

	double C_mean[n_blocks];
	double C_err[n_blocks];
	double P_mean[n_blocks];
	double P_err[n_blocks];

	//campionamento diretto

	double *z_direct = new double[n_throws];

	FillRandomGauss(z_direct, rnd, n_throws, 0, 1);

	for(int i = 0; i < n_throws; i++){

		double S = S0*exp((r - 0.5*sigma*sigma)*T + sigma*z_direct[i]*sqrt(T));

		if(S > K){

			C0[i] = exp(-r*T)*(S - K);
			P0[i] = 0.;

		}else{

			C0[i] = 0.;
			P0[i] = exp(-r*T)*(K - S);

		}
	}

	BlockMethod(C_mean, C_err, C0, n_blocks, L);
	BlockMethod(P_mean, P_err, P0, n_blocks, L);

	delete[] z_direct;

	ofstream fout1("direct.txt");

	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout1 << C_mean[i] << " " << C_err[i] << " " << P_mean[i] << " " << P_err[i] << endl;

	fout1.close();

	//campionamento discreto

	int n_step = 100;
	double dt = T/n_step;

	double *z_discrete = new double[n_throws*n_step];

	FillRandomGauss(z_discrete, rnd, n_throws*n_blocks, 0, 1);

	for(int i = 0; i < n_throws; i++){

		double S = S0;

		for(int j = 0; j < n_step; j++){

			int k = j + i*n_step;
			S *= exp((r - 0.5*sigma*sigma)*dt + sigma*z_discrete[k]*sqrt(dt));

		}

		if(S > K){

			C0[i] = exp(-r*T)*(S - K);
			P0[i] = 0.;

		}else{

			C0[i] = 0.;
			P0[i] = exp(-r*T)*(K - S);

			}
	}

	BlockMethod(C_mean, C_err, C0, n_blocks, L);
	BlockMethod(P_mean, P_err, P0, n_blocks, L);

	delete[] z_discrete;
	delete[] C0;
	delete[] P0;

	ofstream fout2("discrete.txt");

	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout2 << C_mean[i] << " " << C_err[i] << " " << P_mean[i] << " " << P_err[i] << endl;

	fout2.close();

	rnd.SaveSeed(); //salvo il seme della sequenza

	return 0;

}
