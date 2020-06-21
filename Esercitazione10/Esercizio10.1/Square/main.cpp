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

	int l = 32; //numero città
	int temp_step = 100; //numero step MC a temperatura fissata

	double temp = 20; //temperatura di partenza
	double rate = 0.99; //rate di raffreddamento
	double temp_absolute = 0.00001; //temperatura minima

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

	Path *path = new Path();
	path->SetRandom(l,rnd);
	path->SetLenght(mapping);

	cout << "Genero il campione iniziale di traiettorie in modo randomico..." << endl;

	int istep = 0;

	//fase evolutiva
	do{

		double best_lenght = path->GetLenght();

		if(istep%temp_step == 0){

			temp *= rate;
			cout << "Temperatura: " << temp << endl;
			WriteBest(path,temp);		
	
		}

		if(temp < 1.0)

			temp_step = 10000;

		Annealing(l,path,rnd,mapping,temp);
		istep++;
	
	}while(temp > temp_absolute);

	cout << "Fine evoluzione!" << endl;

	WriteFinal(path,l,mapping); //scrivo su file le coordinate del percorso finale

	rnd->SaveSeed(); //salvo il seme della sequenza

	return 0;

}
