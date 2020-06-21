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

	int n_throws = 100000;
	int n_blocks = 100;
	int L = n_throws/n_blocks;
	int n_hit = 0;
	double l = 0.9;
	double d = 1;

	Needle needle[n_throws]; //genero i miei aghi
	
	ThrowNeedle(needle, rnd, d, l, n_throws);

	for(int i = 0; i < n_throws; i++)

		n_hit += needle[i].GetHit();

	double mean[n_blocks];
	double err[n_blocks];

	BlockMethod(mean, err, needle, n_blocks, L, d, l);

	rnd.SaveSeed(); //salvo il seme della sequenza

	ofstream fout("pigreco.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout << mean[i] << " " << err[i] << endl;

	fout.close();

	return 0;

}
