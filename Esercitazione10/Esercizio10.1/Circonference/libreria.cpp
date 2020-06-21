#include "libreria.h"

int PBC(int k, int N){

	if(k >= N)
	
		return k-N;

	else

		return k;

}

double Boltzmann(double temp, double Lold, double Lnew){

	double beta = 1.0/temp;
	double delta = Lnew - Lold;

	return exp(-beta*delta);

}

void Annealing(int l, Path *path, Random *rnd, map<int,City> mapping, double temp){

	double pm_pair = 0.1; //probabilità di scambio coppie
	double pm_shift = 0.1; //probabilità di shift
	double pm_perm = 0.1; //probabilità di permutazioni contigue
	double pm_inversion = 0.1; //probabilità di inversione

	//genero oggetto ausiliario per evoluzione
	Path *pathnew = new Path();
	pathnew->SetRandom(l,rnd);
	pathnew->SetLenght(mapping);
	Equal(pathnew,path,mapping);


	//fase di mutazione genetica

	double p = rnd->Rannyu();

	if(p < pm_pair){

		int i1 = (int)rnd->Rannyu(1,l);
		int i2 = (int)rnd->Rannyu(1,l);

		pathnew->Pair(i1,i2);

	}

	p = rnd->Rannyu();

	if(p < pm_shift){

		int n = (int)rnd->Rannyu(1,l/2);
		int m = (int)rnd->Rannyu(1,l-1);

		pathnew->Shift(n,m);

	}

	p = rnd->Rannyu();

	if(p < pm_perm){

		int m = (int)rnd->Rannyu(2,l/2);

		pathnew->Perm(m);

	}

	p = rnd->Rannyu();

	if(p < pm_inversion){

		int m = (int)rnd->Rannyu(2,l);

		pathnew->Inversion(m);

	}

	pathnew->SetLenght(mapping); //setto la nuova lunghezza per il controllo
	double Lold = path->GetLenght();
	double Lnew = pathnew->GetLenght();

	//accettazione mossa tramite Metropolis
	if(Lnew <= Lold){
	
		Equal(path,pathnew,mapping);

	}else{

		double p = Boltzmann(temp,Lold,Lnew);
		double r = rnd->Rannyu();

		if(r < p)

			Equal(path,pathnew,mapping);

	}

}

void Order(int N, Path *path){

	for(int i = 0; i < N-1; i++){

		int imin = i;
		double min = path[i].GetLenght();

		for(int j = i+1; j < N; j++){

			double l = path[j].GetLenght();

			if(l < min){

				imin = j;
				min = l;
		
			}

		}

		Path appo;
		appo = path[i];
		path[i] = path[imin];
		path[imin] = appo;

	}
}

void Equal(Path *path1, Path *path2, map<int,City> mapping){

	int l = path2->Getl();

	//copio path2 in path1
	for(int i = 0; i < l; i++){

		int c = path2->GetGene(i);
		path1->SetGene(i,c);

	}

	path1->Setl(l);
	path1->SetLenght(mapping);
}

void WriteBest(Path *path, double temp){

	ofstream file;
	file.open("best.path.0",ios::app);
	double best_lenght = path->GetLenght();

	file << temp << " " << best_lenght << endl;
	file.close();

}

void WriteFinal(Path *path, int l, map<int,City> mapping){

	ofstream file;
	file.open("final.path.0");
	map<int,City>::iterator itr;

	for(int i = 0; i < l+1; i++){

		itr = mapping.find(path->GetGene(PBC(i,l)));
		City c = itr->second;

		double x = c.GetX();
		double y = c.GetY();

		file << i << " " << x << " " << y << endl;

	}

	file.close();

}
