#include "libreria.h"

double Error(double *ave, double *ave2, int N){

	if(N == 0)
	
		return 0;

	else
	
		return sqrt((ave2[N] - ave[N]*ave[N])/N);


}

double Psi_100(Point3D x){

	double r = x.GetR();

	return (1.0/pi)*exp(-2.0*r);

}

double Psi_210(Point3D x){

	double r = x.GetR();
	double theta = x.GetTheta();

	return (1.0/(32.0*pi))*r*r*exp(-r)*cos(theta)*cos(theta);

}

void MoveUnif_100(Point3D& point, Random& rnd, int& accepted){

	double x = point.GetX();
	double y = point.GetY();
	double z = point.GetZ();
	double step = point.GetStep();

	double xnew = rnd.Rannyu(x-step,x+step);
	double ynew = rnd.Rannyu(y-step,y+step);
	double znew = rnd.Rannyu(z-step,z+step);
	Point3D pointnew(xnew,ynew,znew);

	double t = Psi_100(pointnew)/Psi_100(point);
	double alpha;
	
	if(t > 1)
		
		alpha = 1;

	else

		alpha = t;

	double s = rnd.Rannyu();

	if(s < alpha){

		point.SetCoord(xnew,ynew,znew);
		accepted++;

	}
}

void MoveUnif_210(Point3D& point, Random& rnd, int& accepted){

	double x = point.GetX();
	double y = point.GetY();
	double z = point.GetZ();
	double step = point.GetStep();

	double xnew = rnd.Rannyu(x-step,x+step);
	double ynew = rnd.Rannyu(y-step,y+step);
	double znew = rnd.Rannyu(z-step,z+step);
	Point3D pointnew(xnew,ynew,znew);

	double t = Psi_210(pointnew)/Psi_210(point);
	double alpha;
	
	if(t > 1)
		
		alpha = 1;

	else

		alpha = t;

	double s = rnd.Rannyu();

	if(s < alpha){

		point.SetCoord(xnew,ynew,znew);
		accepted++;

	}
}

void MoveGauss_100(Point3D& point, Random& rnd, int& accepted){

	double x = point.GetX();
	double y = point.GetY();
	double z = point.GetZ();
	double step = point.GetStep();

	double xnew = rnd.Gauss(x,step);
	double ynew = rnd.Gauss(y,step);
	double znew = rnd.Gauss(z,step);
	Point3D pointnew(xnew,ynew,znew);

	double t = Psi_100(pointnew)/Psi_100(point);
	double alpha;
	
	if(t > 1)
		
		alpha = 1;

	else

		alpha = t;

	double s = rnd.Rannyu();

	if(s < alpha){

		point.SetCoord(xnew,ynew,znew);
		accepted++;

	}
}

void MoveGauss_210(Point3D& point, Random& rnd, int& accepted){

	double x = point.GetX();
	double y = point.GetY();
	double z = point.GetZ();
	double step = point.GetStep();

	double xnew = rnd.Gauss(x,step);
	double ynew = rnd.Gauss(y,step);
	double znew = rnd.Gauss(z,step);
	Point3D pointnew(xnew,ynew,znew);

	double t = Psi_210(pointnew)/Psi_210(point);
	double alpha;
	
	if(t > 1)
		
		alpha = 1;

	else

		alpha = t;

	double s = rnd.Rannyu();

	if(s < alpha){

		point.SetCoord(xnew,ynew,znew);
		accepted++;

	}
}

void EquilUnif_100(int nstep, Point3D point, double r_real, Random rnd){

	bool restart = true;
	int count = 0; //tengo conto del numero di cicli necessari all'equilibrio
	int accepted = 0;

	do{

		double sum = 0., sum2 = 0.;
		double r, err;

		for(int i = 0; i < nstep; i++){

			MoveUnif_100(point,rnd,accepted);
			sum += point.GetR();
			sum2 += point.GetR()*point.GetR();

		}

		sum /= nstep;
		sum2 /= nstep;
		r = sum;
		err = sqrt((sum2 - sum*sum));

		if(abs(r-r_real) < err) restart = false;

		count ++;

	}while(restart == true);

	cout << "Funzione d'onda stato fondamentale con probabilità di transizione uniforme" << endl;
	cout << "Per raggiungere la condizione d'equilibrio sono stati necessari: " << count << endl << endl;

}

void EquilUnif_210(int nstep, Point3D point, double r_real, Random rnd){

	bool restart = true;
	int count = 0; //tengo conto del numero di cicli necessari all'equilibrio
	int accepted = 0;

	do{

		double sum = 0., sum2 = 0.;
		double r, err;

		for(int i = 0; i < nstep; i++){

			MoveUnif_210(point,rnd,accepted);
			sum += point.GetR();
			sum2 += point.GetR()*point.GetR();

		}

		sum /= nstep;
		sum2 /= nstep;
		r = sum;
		err = sqrt((sum2 - sum*sum));

		if(abs(r-r_real) < err) restart = false;

		count ++;

	}while(restart == true);

	cout << "Funzione d'onda stato eccitato con probabilità di transizione uniforme" << endl;
	cout << "Per raggiungere la condizione d'equilibrio sono stati necessari: " << count << endl << endl;

}

void EquilGauss_100(int nstep, Point3D point, double r_real, Random rnd){

	bool restart = true;
	int count = 0; //tengo conto del numero di cicli necessari all'equilibrio
	int accepted = 0;

	do{

		double sum = 0., sum2 = 0.;
		double r, err;

		for(int i = 0; i < nstep; i++){

			MoveGauss_100(point,rnd,accepted);
			sum += point.GetR();
			sum2 += point.GetR()*point.GetR();

		}

		sum /= nstep;
		sum2 /= nstep;
		r = sum;
		err = sqrt((sum2 - sum*sum));

		if(abs(r-r_real) < err) restart = false;

		count ++;

	}while(restart == true);

	cout << "Funzione d'onda stato fondamentale con probabilità di transizione normale (gaussiana)" << endl;
	cout << "Per raggiungere la condizione d'equilibrio sono stati necessari: " << count << endl << endl;

}

void EquilGauss_210(int nstep, Point3D point, double r_real, Random rnd){

	bool restart = true;
	int count = 0; //tengo conto del numero di cicli necessari all'equilibrio
	int accepted = 0;

	do{

		double sum = 0., sum2 = 0.;
		double r, err;

		for(int i = 0; i < nstep; i++){

			MoveGauss_210(point,rnd,accepted);
			sum += point.GetR();
			sum2 += point.GetR()*point.GetR();

		}

		sum /= nstep;
		sum2 /= nstep;
		r = sum;
		err = sqrt((sum2 - sum*sum));

		if(abs(r-r_real) < err) restart = false;

		count ++;

	}while(restart == true);

	cout << "Funzione d'onda stato eccitato con probabilità di transizione normale (gaussiana)" << endl;
	cout << "Per raggiungere la condizione d'equilibrio sono stati necessari: " << count << endl << endl;

}

void DataUnif_100(int n_throws, Point3D point, double *r, Random rnd){

	ofstream fout("orbital100.txt");
	int accepted = 0;

	for(int i = 0; i < n_throws; i++){

		MoveUnif_100(point,rnd,accepted);
		r[i] = point.GetR();

		if(i%100 == 0)

			fout << point.GetX() << " " << point.GetY() << " " << point.GetZ() << endl;
	}

	fout.close();
	cout << "Rate di accettazione = " << (double)accepted/n_throws << endl << endl;
}

void DataUnif_210(int n_throws, Point3D point, double *r, Random rnd){

	ofstream fout("orbital210.txt");
	int accepted = 0;

	for(int i = 0; i < n_throws; i++){

		MoveUnif_210(point,rnd,accepted);
		r[i] = point.GetR();

		if(i%100 == 0)

			fout << point.GetX() << " " << point.GetY() << " " << point.GetZ() << endl;	
	}

	fout.close();
	cout << "Rate di accettazione = " << (double)accepted/n_throws << endl << endl;
}

void DataGauss_100(int n_throws, Point3D point, double *r, Random rnd){

	int accepted = 0;

	for(int i = 0; i < n_throws; i++){

		MoveGauss_100(point,rnd,accepted);
		r[i] = point.GetR();	

	}
	cout << "Rate di accettazione = " << (double)accepted/n_throws << endl << endl;
}

void DataGauss_210(int n_throws, Point3D point, double *r, Random rnd){

	int accepted = 0;

	for(int i = 0; i < n_throws; i++){

		MoveGauss_210(point,rnd,accepted);
		r[i] = point.GetR();	

	}
	cout << "Rate di accettazione = " << (double)accepted/n_throws << endl << endl;
}

void BlockMethod(double *mean, double *err_mean, double *data, int n_blocks, int L){

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
			accu += data[k];

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
