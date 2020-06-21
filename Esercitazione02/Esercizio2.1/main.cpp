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

	Function1D * f = new Cos1D();

	int n_throws = 10000;
	int n_blocks = 100;
	int L = n_throws/n_blocks;

	//prima parte: campiono una distribuzione uniforme nell'intervallo (0,1)

	double *x1 = new double[n_throws*n_blocks]; //vettore contenente i numeri casuali necessari all'esperimento
						    //alloco dinamicamente in modo da avere a disposizione più memoria

	double mean_UNIF[n_blocks];
	double err_UNIF[n_blocks];

	UniformSampling(x1, rnd, n_throws*n_blocks);

	BlockMethod(mean_UNIF, err_UNIF, x1, f, n_blocks, L);

	ofstream fout1("uniform.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout1 << mean_UNIF[i] << " " << err_UNIF[i] << endl;

	fout1.close();

	delete [] x1;

	//seconda parte: (importance sampling) campiono una distribuzione non-uniforme nell'intervallo (0,1)
	//come distribuzione di probabilità, utilizzo una retta passante per i punti (0,pi/2) e (1,0), campionandola nell'intervallo (0,1)
	//utilizzandando una distribuzione semplice come la retta, è possibile applicare il metodo d'inversione della cumulativa

	Function1D * g = new ImportanceFun(); //nuova funzione ottenuta con il metodo importance sampling

	double *x2 = new double[n_throws*n_blocks]; //vettore contenente i numeri casuali necessari all'esperimento
						    //alloco dinamicamente in modo da avere a disposizione più memoria

	double mean_IMP[n_blocks];
	double err_IMP[n_blocks];

	ImportanceSampling(x2, rnd, n_throws*n_blocks); //riempio il vettore di numeri casuali

	BlockMethod(mean_IMP, err_IMP, x2, g, n_blocks, L);

	ofstream fout2("importance.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout2 << mean_IMP[i] << " " << err_IMP[i] << endl;

	fout2.close();

	delete [] x2;

	rnd.SaveSeed(); //salvo il seme della sequenza

	return 0;

}
