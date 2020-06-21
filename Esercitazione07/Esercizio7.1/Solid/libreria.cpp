#include "libreria.h"

double Mean(double *x, int n, int t){

	double sum = 0.0, m = 0.0;

	for(int i = 0; i < n-t; i+=t){

		sum += x[i];
		m += 1;

	}

	return sum/(double)m;


}

double Variance(double *x, int n, int t){

	double sum = 0.0, sum2 = 0.0, m = 0.0;

	for(int i = 0; i < n-t; i+=t){

		sum += x[i];
		sum2 += x[i]*x[i];
		m += 1;

	}

	sum /= (double)m;
	sum2 /= (double)m;

	return sum2 - sum*sum;


}

double Autocorr(double *x, int n, int t){

	double s = 0.0, m = 0.0;

	double mean = Mean(x,n,t);
	double var = Variance(x,n,t);

	for(int i = t; i < n; i += t){

		s += x[i-t]*x[i];
		m += 1;

	}

	s /= m;
	
	return (s - mean*mean)/var;

}

double Error(double ave, double ave2, int n){

	if(n == 0)
	
		return 0;

	else
	
		return sqrt((ave2 - ave*ave)/(n-1));


}

double ErrorBM(double *x, int N, int L){

	int nblk = N/L;

	//genero vettori ausiliari
	double ave[nblk];
	double ave2[nblk];
	double sum = 0.0;
	double sum2 = 0.0;

	for(int i = 0; i < nblk; i++){

		//azzero le componenti dei vettori di accumulazione
		ave[i] = 0;
		ave2[i] = 0;

	}


	for(int i = 0; i < nblk; i++){

		double accu = 0;

		for(int j = 0; j < L; j++){

			int k = j + i*L;
			accu += x[k];

		}

		ave[i] = accu/L; //valore d'aspettazione per ogni blocco
		ave2[i] = ave[i]*ave[i]; //varianza associata ad ogni blocco

	}

	for(int j = 0; j < nblk; j++){

		sum += ave[j];
		sum2 += ave2[j];

	}

	sum /= nblk;
	sum2 /= nblk;

	return Error(sum,sum2,nblk);

}
