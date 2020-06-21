#include "libreria.h"

int main(int argc, char *argv[]){

	int N = 500000; //numero dati
	double *epot, *pres;

	epot = new double[N];
	pres = new double[N];

	//carico da file i dati di energia interna e pressione istantanei
	ifstream Epot, Pres;

	Epot.open("Data/ist_epot.out");
	Pres.open("Data/ist_pres.out");

	for(int i = 0; i < N; i++){

		Epot >> epot[i];
		Pres >> pres[i];

	}

	Epot.close();
	Pres.close();

	ofstream acorr_epot, acorr_pres;

	acorr_epot.open("acorr_epot.out");
	acorr_pres.open("acorr_pres.out");

	int M = 100;

	for(int i = 0; i < M; i++){

		double U_acorr = Autocorr(epot, N, i+1);
		double P_acorr = Autocorr(pres, N, i+1);

		acorr_epot << i+1 << " " << U_acorr << endl;
		acorr_pres << i+1 << " " << P_acorr << endl;

	}

	acorr_epot.close();
	acorr_pres.close();

	M = 500;

	ofstream err_epot, err_pres;

	err_epot.open("err_epot.out");
	err_pres.open("err_pres.out");

	for(int i = 0; i < M; i++){

		double U_err = ErrorBM(epot, N, 10*(i+1));
		double P_err = ErrorBM(pres, N, 10*(i+1));

		err_epot << 10*(i+1) << " " << U_err << endl;
		err_pres << 10*(i+1) << " " << P_err << endl;

	}

	err_epot.close();
	err_pres.close();

	return 0;

}
