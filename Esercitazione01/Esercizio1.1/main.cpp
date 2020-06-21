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

	int M = 100000; //numero totale di lanci
	int N = 100; //numero di blocchi (Block Method)
	int L = int(M/N); //numero lanci per ogni blocco

	//Prima Parte dell'esercizio

	double r1[M]; //vettore contenente i numeri casuali necessari all'esperimento

	double mean[N];
	double err_mean[N];

	//carico il vettore r con numeri (pseudo)casuali distribuiti uniformemente tra 0 e 1
	FillRandom(r1, rnd, M, 0, 1);
	rnd.SaveSeed(); //salvo il seme della sequenza, per garantire la riproducibilità della sequenza

	BlockMean(mean, err_mean, r1, N, L);

	ofstream fout1("BM-mean.txt");

	
	for (int i = 0; i < N; i++)

		//carico su file i risultati
		fout1 << mean[i] << " " << err_mean[i] << endl;

	fout1.close();

	//Seconda Parte dell'esercizio

	double var[N];
	double err_var[N];

	BlockVariance(var, err_var, r1, N, L);

	ofstream fout2("BM-variance.txt");
	
	for (int i = 0; i < N; i++)

		//carico su file i risultati
		fout2 << var[i] << " " << err_var[i] << endl;

	fout2.close();
		
	//Terza Parte: Test Chi Quadro

	int n_step = 100; //numero sottointervalli nell'intervallo (0,1)
	int n_throws = 10000; //numero lanci per blocco
	double *r2 = new double[n_throws*n_step]; //vettore contenente i numeri casuali necessari all'esperimento
						  //alloco dinamicamente in modo da avere a disposizione più memoria

	FillRandom(r2, rnd, n_throws*n_step, 0, 1); //genero 10^6 numeri (pseudo)casuali

	int n[n_step]; //vettore contenente il numero di conteggi in ogni sottointervallo

	ofstream fout3("chiquadro.txt");

	for(int i = 0; i < n_step; i++){

		double chi = ChiQuadro(n, r2, n_step, n_throws, i);
		//carico su file i risultati
		fout3 << chi << endl;

	}

	delete [] r2;

	fout3.close();

	return 0;

}
