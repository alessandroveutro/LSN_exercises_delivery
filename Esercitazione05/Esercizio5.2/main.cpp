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

	int n_throws = 1000000; //numero step
	int n_blocks = 100; //numero blocchi
	int nstep_eq = 100;
	int L = (int)n_throws/n_blocks;

	double step;

	double *r;
	r = new double[n_throws];

	double r100_mean[n_blocks], r100_err[n_blocks];
	double r210_mean[n_blocks], r210_err[n_blocks];

	double r100_real = 1.5;
	double r210_real = 5.0;

	Point3D point(0.,0.,0.); //setto la classe punto 3D

	//Uniform Transition Probability
	//Funzione d'onda stato fondamentale

	step = 0.3; //setto inizialmente uno step casuale
	point.SetStep(step); //setto lo step	
	point.SetCoord(1.5,0.,0.); //setto il punto iniziale
	EquilUnif_100(nstep_eq, point, r100_real, rnd); //fase di equilibrazione

	step = 1.2; //setto uno step scelto empiricamente per ottenere un'accettazione del 50%
	point.SetStep(step); //setto lo step
	DataUnif_100(n_throws, point, r, rnd);
	BlockMethod(r100_mean, r100_err, r, n_blocks, L);

	ofstream fout1("unif100.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout1 << r100_mean[i] << " " << r100_err[i] << endl;

	fout1.close();

	//Funzione d'onda stato eccitato

	step = 0.5; //setto inizialmente uno step casuale
	point.SetStep(step); //setto lo step	
	point.SetCoord(5.,0.,0.); //setto il punto iniziale
	EquilUnif_210(nstep_eq, point, r210_real, rnd); //fase di equilibrazione

	step = 3.; //setto uno step scelto empiricamente per ottenere un'accettazione del 50%
	point.SetStep(step); //setto lo step
	DataUnif_210(n_throws, point, r, rnd);
	BlockMethod(r210_mean, r210_err, r, n_blocks, L);

	ofstream fout2("unif210.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout2 << r210_mean[i] << " " << r210_err[i] << endl;

	fout2.close();

	//Normal Transition Probability
	//Funzione d'onda stato fondamentale

	step = 0.2; //setto inizialmente uno step casuale
	point.SetStep(step); //setto lo step	
	point.SetCoord(10,0.,0.); //setto il punto iniziale
	EquilGauss_100(nstep_eq, point, r100_real, rnd); //fase di equilibrazione

	step = 0.7; //setto uno step scelto empiricamente per ottenere un'accettazione del 50%
	point.SetStep(step); //setto lo step
	DataGauss_100(n_throws, point, r, rnd);
	BlockMethod(r100_mean, r100_err, r, n_blocks, L);

	ofstream fout3("gauss100.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout3 << r100_mean[i] << " " << r100_err[i] << endl;

	fout3.close();

	//Funzione d'onda stato eccitato

	step = 0.2; //setto inizialmente uno step casuale
	point.SetStep(step); //setto lo step	
	point.SetCoord(5.,0.,0.); //setto il punto iniziale
	EquilGauss_210(nstep_eq, point, r210_real, rnd); //fase di equilibrazione

	step = 1.7; //setto uno step scelto empiricamente per ottenere un'accettazione del 50%
	point.SetStep(step); //setto lo step
	DataGauss_210(n_throws, point, r, rnd);
	BlockMethod(r210_mean, r210_err, r, n_blocks, L);

	ofstream fout4("gauss210.txt");
	
	for (int i = 0; i < n_blocks; i++)

		//carico su file i risultati
		fout4 << r210_mean[i] << " " << r210_err[i] << endl;

	fout4.close();

	delete[] r;

	rnd.SaveSeed(); //salvo il seme della sequenza

	return 0;

}
