#include "libreria.h"


double Error(double *ave, double *ave2, int N){

	if(N == 0)
	
		return 0;

	else
	
		return sqrt((ave2[N] - ave[N]*ave[N])/N);


}


void ThrowNeedle(Needle *needle, Random rnd, double d, double l, int n_throws){

	for(int i = 0; i < n_throws; i ++){

		//Genero la coordinata y del centro del mio ago sul piano (x,y) all'interno di una griglia (10*d, 10*d)
		double yc = rnd.Rannyu(d, 10*d);

		//per ottenere l'angolo di inclinazione dell'ago distribuito in modo uniforme in (0,2*pi), utilizzo la tecnica del rigetto:
		//genero punti all'interno di una circonferenza di raggio unitario e inverto la legge di trasformazione delle coordinate sferiche,
		//in modo da avere 'theta = theta(x,y)'
		double theta;

		for(int j = 0; j < 1; j++){

			double x = rnd.Rannyu();
			double y = rnd.Rannyu();
			double r = sqrt(x*x + y*y);

			if(r < 1){

				if(y >= 0){
		
					theta = acos(x/r);

				}else{

					theta = -acos(x/r);

				}

			}else
		
				j--;

		}

		needle[i].Throw(yc,l,theta);
		needle[i].CheckHit(d);		

	}
}


void BlockMethod(double *mean, double *err_mean, Needle *needle, int n_blocks, int L, double d, double l){

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

		int n_hit = 0;

		for(int j = 0; j < L; j++){

			int k = j + i*L;
			n_hit += needle[k].GetHit(); //sommo il numero di successi per ogni blocco

		}

		ave[i] = 2*L*l/(n_hit*d);
		ave2[i] = ave[i]*ave[i];

	}

	for(int i = 0; i < n_blocks; i++){

		for(int j = 0; j < i + 1; j++){

			sum[i] += ave[j];
			sum2[i] += ave2[j];

		}

		sum[i] /= (i + 1);
		sum2[i] /= (i + 1);
		mean[i] = sum[i];
		err_mean[i] = Error(sum, sum2, i); //calcolo l'errore su i+1 misure

	}
}
