#include "libreria.h"


double Error(double *ave, double *ave2, int N){

	if(N == 0)
	
		return 0;

	else
	
		return sqrt((ave2[N] - ave[N]*ave[N])/N);


}

double N(double x){

	return 0.5*(1 + erf(x/sqrt(2)));

}

double C_BlackScholes(double S0, double K, double T, double r, double sigma){

	double d1 = (1/(sigma*sqrt(T)))*(log(S0/K) + (r + 0.5*sigma*sigma*T));
	double d2 = d1 - sigma*sqrt(T);

	return S0*N(d1) - K*exp(-r*T)*N(d2);

}

double P_BlackScholes(double S0, double K, double T, double r, double sigma){

	double d1 = (1/(sigma*sqrt(T)))*(log(S0/K) + (r + 0.5*sigma*sigma*T));
	double d2 = d1 - sigma*sqrt(T);

	return S0*(N(d1)-1) - K*exp(-r*T)*(N(d2)-1);

}

void FillRandomGauss(double *z, Random rnd, int dim, double mean, double sigma){

	for(int i = 0; i < dim; i++)

		z[i] = rnd.Gauss(mean,sigma);

}

void BlockMethod(double *mean, double *err_mean, double *x, int n_blocks, int L){

	//genero vettori ausiliari
	double ave[n_blocks];
	double ave2[n_blocks];
	double sum[n_blocks];
	double sum2[n_blocks];

	for(int i = 0; i < n_blocks; i++){

		//azzero le componenti dei vettori di accumulazione
		ave[i] = 0;
		ave2[i] = 0;
		sum[i] = 0;
		sum2[i] = 0;

	}


	for(int i = 0; i < n_blocks; i++){

		double accu = 0;

		for(int j = 0; j < L; j++){

			int k = j + i*L;
			accu += x[k];

		}

		ave[i] = accu/L; //valore d'aspettazione per ogni blocco
		ave2[i] = ave[i]*ave[i]; //varianza associata ad ogni blocco

	}

	for(int i = 0; i < n_blocks; i++){

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
