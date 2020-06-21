#include "libreria.h"


double Error(double *ave, double *ave2, int N){

	if(N == 0)
	
		return 0;

	else
	
		return sqrt((ave2[N] - ave[N]*ave[N])/N);


}

void FillRandom(double *r, Random rnd, int dim, double min, double max){

	for(int i = 0; i < dim; i++)

		r[i] = rnd.Rannyu(min, max);

}

void BlockMean(double *mean, double *err_mean, double *r, int N, int L){

	//genero vettori ausiliari
	double ave[N];
	double ave2[N];
	double sum[N];
	double sum2[N];

	for(int i = 0; i < N; i++){

		//azzero le componenti dei vettori di accumulazione
		ave[i] = 0;
		ave2[i] = 0;
		sum[i] = 0;
		sum2[i] = 0;

	}


	for(int i = 0; i < N; i++){

		double accu = 0;

		for(int j = 0; j < L; j++){

			int k = j + i*L;
			accu += r[k];

		}

		ave[i] = accu/L; //valore d'aspettazione per ogni blocco
		ave2[i] = ave[i]*ave[i]; //varianza associata ad ogni blocco

	}

	for(int i = 0; i < N; i++){

		for(int j = 0; j < i + 1; j++){

			sum[i] += ave[j];
			sum2[i] += ave2[j];

		}

		sum[i] /= (i + 1);
		sum2[i] /= (i + 1);
		mean[i] = sum[i]; //valore medio al passo i+1
		err_mean[i] = Error(sum, sum2, i); //incertezza sul valore medio al passo i+1

	}


}

void BlockVariance(double *var, double *err_var, double *r, int N, int L){

	//genero vettori ausiliari
	double ave[N];
	double ave2[N];
	double sum[N];
	double sum2[N];

	for(int i = 0; i < N; i++){

		//azzero le componenti dei vettori di accumulazione
		ave[i] = 0;
		ave2[i] = 0;
		sum[i] = 0;
		sum2[i] = 0;

	}


	for(int i = 0; i < N; i++){

		double accu = 0;

		for(int j = 0; j < L; j++){

			int k = j + i*L;
			accu += (r[k] - 0.5)*(r[k] - 0.5);

		}

		ave[i] = accu/L; //valore d'aspettazione per ogni blocco
		ave2[i] = ave[i]*ave[i]; //varianza associata ad ogni blocco

	}

	for(int i = 0; i < N; i++){

		for(int j = 0; j < i + 1; j++){

			sum[i] += ave[j];
			sum2[i] += ave2[j];

		}

		sum[i] /= (i + 1);
		sum2[i] /= (i + 1);
		var[i] = sum[i]; //valore medio della deviazione standard al passo i+1
		err_var[i] = Error(sum, sum2, i); //calcolo l'incertezza sulla deviazione standard al passo i+1

	}
}

double ChiQuadro(int *n, double *r, int n_step, int n_throws, int l){

	for(int i = 0; i < n_step; i++)

		n[i] = 0; //azzero le componenti

	for(int i = 0; i < n_step; i++)

		for(int j = 0; j < n_throws; j++){

			int k = j + l*n_throws;
			if(r[k] >= (double)i/n_step and r[k] < (double)(i+1)/n_step)

				n[i]++;

		}

	double accu = 0;
	double mean = (double)n_throws/n_step;
	
	for(int i = 0; i < n_step; i++)

		accu += (n[i] - mean)*(n[i] - mean);

	return accu/mean;

}
