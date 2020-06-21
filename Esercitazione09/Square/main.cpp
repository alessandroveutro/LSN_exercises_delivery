#include "libreria.h"

#define Pi 3.14152629

int main(int argc, char *argv[]){

	//setto il generatore di numeri casuali
	Random *rnd;
	rnd = new Random(); 
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
			rnd->SetRandom(seed,p1,p2);
		}
	}

	input.close();

	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	//inizio esercizio

	int N = 100; //popolazione
	int l = 32; //numero città
	int nsteps = 100000; //numero step di simulazione

	double a = 1.0;

	map<int,City> mapping;

	//genero N città all'interno di un quadrato
	for(int i = 1; i < l+1; i++){

		double x = rnd->Rannyu(0,a);
		double y = rnd->Rannyu(0,a);
		City c;
		c.Set(x,y);
		mapping.insert(pair<int,City>(i,c));

	}

	Path path[N];

	cout << "Simulazione algoritmo genetico" << endl;
	cout << "Numero città da visitare: " << l << endl;
	cout << "Popolazione campione: " << N << endl;
	cout << "Step evolutivi: " << nsteps << endl;

	cout << "Genero il campione iniziale di traiettorie in modo randomico..." << endl;
	//genero il campione iniziale di traiettorie in modo randomico
	for(int i = 0; i < N; i++){

		path[i].SetRandom(l,rnd);
		path[i].SetLenght(mapping); //setto la lunghezza di ogni traiettoria

	}

	Order(N,path); //Ordino il mio campione per lunghezza della traiettoria

	//fase evolutiva
	for(int i = 0; i < nsteps; i++){

		double sum = 0.0, sum2 = 0.0, error;
		double best_lenght = path[0].GetLenght();

		for(int j = 0; j < 50; j++){

			double lenght = path[j].GetLenght();
			sum += lenght;
			sum2 += lenght*lenght;

		}

		sum /= (double)50;
		sum2 /= (double)50;
		error = sqrt(abs(sum2-sum*sum));

		WriteBest(path[0],i);
		WriteMean(path,i);
		Evolution(N,l,path,rnd,mapping);

	}

	cout << "Fine evoluzione!" << endl;
	Order(N,path); //Ordino il mio campione per lunghezza della traiettoria

	WriteFinal(path[0],l,mapping); //scrivo su file le coordinate del percorso finale

	map<int,City>::iterator itr;

	for(int i = 0; i < l+1; i++){

		itr = mapping.find(path[0].GetGene(PBC(i,l)));
		City c = itr->second;

	}

	rnd->SaveSeed(); //salvo il seme della sequenza

	return 0;

}
