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

	int n_walks = 100000; //numero random walk
	int n_steps = 100; //numero passi
	double a = 1; //passo del reticolo

	RandomWalk RW[n_walks]; //classe oggetto che descrive un generico random walk

	//prima parte: random walk discrete

	double *random1 = new double[n_steps*n_walks];
	FillRandomDiscrete(random1, rnd, n_steps, n_walks);

	double mean_disc[n_steps];
	double err_disc[n_steps];

	RW_Discrete(mean_disc, err_disc, RW, random1, a, n_steps, n_walks);

	ofstream fout1("discrete.txt");
	
	for (int i = 0; i < n_steps; i++)

		//carico su file i risultati
		fout1 << mean_disc[i] << " " << err_disc[i] << endl;

	fout1.close();

	delete[] random1;


	//resetto il random walk
	for (int i = 0; i < n_walks; i++)

		RW[i].Reset();


	//seconda parte: random walk continuum

	double *random2 = new double[n_steps*n_walks];
	double *random3 = new double[n_steps*n_walks];
	FillRandomContinuum(random2, random3, rnd, n_steps, n_walks);

	double mean_cont[n_steps];
	double err_cont[n_steps];

	RW_Continuum(mean_cont, err_cont, RW, random2, random3, a, n_steps, n_walks);

	ofstream fout2("continuum.txt");
	
	for (int i = 0; i < n_steps; i++)

		//carico su file i risultati
		fout2 << mean_cont[i] << " " << err_cont[i] << endl;

	fout2.close();

	delete[] random2;
	delete[] random3;

	rnd.SaveSeed(); //salvo il seme della sequenza

	return 0;

}
