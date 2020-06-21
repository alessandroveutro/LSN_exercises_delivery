#include "random.h"
#include "libreria.h"

int main(int argc, char *argv[]){

	int n_throws = 10000; //numero lanci totali per ogni dado

	Random rnd; //setto il generatore di numeri casuali
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


	double S_2[n_throws];
	double S_10[n_throws];
	double S_100[n_throws];

	double *r = new double[n_throws*100]; //vettore contenente i numeri casuali necessari all'esperimento

	//Standard Dice

	StandardDice(r, rnd, n_throws*100);

	//N = 2
	for(int i = 0; i < n_throws; i++){

		S_2[i] = 0;

		for(int j = 0; j < 2; j++){
	
			int k = i + j*n_throws;
			S_2[i] += r[k];

		}
	
		S_2[i] /= 2;

	}


	//N = 10
	for(int i = 0; i < n_throws; i++){

		S_10[i] = 0;

		for(int j = 0; j < 10; j++){
	
			int k = i + j*n_throws;
			S_10[i] += r[k];

		}
	
		S_10[i] /= 10;
	}


	//N = 100
	for(int i = 0; i < n_throws; i++){

		S_100[i] = 0;

		for(int j = 0; j < 100; j++){
	
			int k = i + j*n_throws;
			S_100[i] += r[k];

		}
	
		S_100[i] /= 100;

	}

	ofstream fout1("standard.txt");

	for(int i = 0; i < n_throws; i++)

		//carico su file i risultati
		fout1 << r[i] << " " << S_2[i] << " " << S_10[i] << " " << S_100[i] << endl;

	fout1.close();

	//Exponential Dice

	double lambda = 1;
	ExpDice(r, rnd, n_throws*100, lambda);

	//N = 2
	for(int i = 0; i < n_throws; i++){

		S_2[i] = 0;

		for(int j = 0; j < 2; j++){
	
			int k = i + j*n_throws;
			S_2[i] += r[k];

		}
	
		S_2[i] /= 2;

	}


	//N = 10
	for(int i = 0; i < n_throws; i++){

		S_10[i] = 0;

		for(int j = 0; j < 10; j++){
	
			int k = i + j*n_throws;
			S_10[i] += r[k];

		}
	
		S_10[i] /= 10;

	}


	//N = 100
	for(int i = 0; i < n_throws; i++){

		S_100[i] = 0;

		for(int j = 0; j < 100; j++){
	
			int k = i + j*n_throws;
			S_100[i] += r[k];

		}
	
		S_100[i] /= 100;

	}


	ofstream fout2("exponential.txt");

	for(int i = 0; i < n_throws; i++)

		//carico su file i risultati
		fout2 << r[i] << " " << S_2[i] << " " << S_10[i] << " " << S_100[i] << endl;

	fout2.close();

	//Cauchy-Lorentz Dice

	double mu = 0;
	double gamma = 1;
	CauchyDice(r, rnd, n_throws*100, mu, gamma);

	//N = 2
	for(int i = 0; i < n_throws; i++){

		S_2[i] = 0;

		for(int j = 0; j < 2; j++){
	
			int k = i + j*n_throws;
			S_2[i] += r[k];

		}
	
		S_2[i] /= 2;

	}


	//N = 10
	for(int i = 0; i < n_throws; i++){

		S_10[i] = 0;

		for(int j = 0; j < 10; j++){
	
			int k = i + j*n_throws;
			S_10[i] += r[k];

		}
	
		S_10[i] /= 10;

	}


	//N = 100
	for(int i = 0; i < n_throws; i++){

		S_100[i] = 0;

		for(int j = 0; j < 100; j++){
	
			int k = i + j*n_throws;
			S_100[i] += r[k];

		}
	
		S_100[i] /= 100;

	}

	ofstream fout3("cauchy.txt");

	for(int i = 0; i < n_throws; i++)

		//carico su file i risultati
		fout3 << r[i] << " " << S_2[i] << " " << S_10[i] << " " << S_100[i] << endl;

	fout3.close();

	delete[] r;

	rnd.SaveSeed(); //salvo il seme della sequenza

	return 0;

}
