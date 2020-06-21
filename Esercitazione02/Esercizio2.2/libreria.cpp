#include "libreria.h"


double Error(double *ave, double *ave2, int N){

	if(N == 0)
	
		return 0;

	else
	
		return sqrt((ave2[N] - ave[N]*ave[N])/N);


}

void FillRandomDiscrete(double *r, Random random, int n_steps, int n_walks){

	int dim = n_steps*n_walks;
	
	for(int i = 0; i < dim; i++)
	
		r[i] = (int)random.Rannyu(0,6);

}


void FillRandomContinuum(double *r1, double *r2, Random random, int n_steps, int n_walks){

	int dim = n_steps*n_walks;
	
	for(int i = 0; i < dim; i ++){
	
		r1[i] = random.Rannyu(0,2*pi);
		r2[i] = random.Rannyu();

	}

}

void RW_Discrete(double *mean, double *err_mean, RandomWalk *rw, double *rand, double a, int n_steps, int n_walks){

	//genero vettori ausiliari
	double sum[n_steps];
	double sum2[n_steps];

	for(int i = 0; i < n_steps; i++){

		//azzero le componenti dei vettori di accumulazione
		sum[i] = 0;
		sum2[i] = 0;

	}


	for(int i = 0; i < n_steps; i++){

		for(int j = 0; j < n_walks; j++){

			int k = j + i*n_walks;
			rw[j].StepDiscrete(rand[k],a); //genero un random walk discreto all'interno di un reticolo con passo a
			double r2 = rw[j].GetR2();
			sum[i] += r2;
			sum2[i] += r2*r2;

		}

		sum[i] /= n_walks; //valor medio del cammino quadratico medio al passo i-esimo
		sum2[i] /= n_walks;
		mean[i] = sqrt(sum[i]); //calcolo la radice del valor medio del cammino quadratico medio al passo i-esimo
		err_mean[i] = 0.5*Error(sum, sum2, i)/mean[i]; //incertezza sul valor medio, calcolata propagando le incertezze

	}
}

void RW_Continuum(double *mean, double *err_mean, RandomWalk *rw, double *rand1, double *rand2, double a, int n_steps, int n_walks){

	//genero vettori ausiliari
	double sum[n_steps];
	double sum2[n_steps];

	for(int i = 0; i < n_steps; i++){

		//azzero le componenti dei vettori di accumulazione
		sum[i] = 0;
		sum2[i] = 0;

	}

	for(int i = 0; i < n_steps; i++){

		for(int j = 0; j < n_walks; j++){

			int k = j + i*n_walks;
			rw[j].StepContinuum(rand1[k],rand2[k],a); //genero un random walk discreto all'interno di un reticolo con passo a
			double r2 = rw[j].GetR2();
			sum[i] += r2;
			sum2[i] += r2*r2;

		}

		sum[i] /= n_walks; //valor medio del cammino quadratico medio al passo i-esimo
		sum2[i] /= n_walks;
		mean[i] = sqrt(sum[i]); //calcolo la radice del valor medio del cammino quadratico medio al passo i-esimo
		err_mean[i] = 0.5*Error(sum, sum2, i)/mean[i]; //incertezza sul valor medio, calcolata propagando le incertezze

	}
}
