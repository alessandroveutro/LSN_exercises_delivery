#include "libreria.h"

int Select(int N, Random *rnd){

	double t = rnd->Rannyu();
	double p = 5.;

	return (int)N*pow(t,p);
}


int PBC(int k, int N){

	if(k >= N)
	
		return k-N;

	else

		return k;

}


void Crossover(Path *path, int i, int N, int l, int n){

	int i1 = i;
	int i2 = i + N/2;

	//vettori contenenti i geni finali della catena
	int n1[n];
	int n2[n];
	int k = 0;

	for(int i = l-n; i < l; i++){

		n1[k] = path[i1].GetGene(i);
		n2[k] = path[i2].GetGene(i);
		k++;

	}

	//fase di incrocio dei campioni
	Cross(path,n1,i1,i2,l,n);
	Cross(path,n2,i2,i1,l,n);

}

void Cross(Path *path, int *genes, int i1, int i2, int l, int n){ //codice genetico i1 viene modificato con quello di i2, che rimane invece invariato

	int ind = l-n; //primo indice da riempire

	for(int i = 0; i < l-n; i++){

		int count = 0;
		int c = path[i2].GetGene(i);

		for(int j = 0; j < l-n; j++){

			if(c != path[i1].GetGene(j))

				count++;

		}

		if(count == l-n){

			path[i1].SetGene(ind,c);
			ind++;

		}
	}

	//se sono ancora riuscito a riempire la mia catena,
	//inserisco il codice genetico mancante seguendo l'ordinamento iniziale
	if(ind != l){

		for(int i = 0; i < n; i++){

			int count = 0;
			int c = genes[i];

			for(int j = l-n; j < ind; j++){

				if(c != path[i1].GetGene(j))

					count++;

			}

			if(count == ind-l+n){

				path[i1].SetGene(ind,c);
				ind++;

			}
		}
	}
}

void Evolution(int N, int l, Path *path, Random *rnd, map<int,City> mapping){

	double pm_pair = 0.1; //probabilità di scambio coppie
	double pm_shift = 0.07; //probabilità di shift
	double pm_perm = 0.1; //probabilità di permutazioni contigue
	double pm_inversion = 0.1; //probabilità di inversione
	double pc = 0.7; //probabilità di cross

	Path pathnew[N];

	//inizializzo il vettore
	for(int i = 0; i < N; i++){

		pathnew[i].SetRandom(l,rnd);

	}

	Order(N,path); //ordino la catena per il crossover

	//fase di crossover genetico
	for(int i = 0; i < N/2; i++){

		//seleziono due campioni 'favoriti' distinti per l'incrocio genetico
		int i1 = Select(N, rnd);
		int i2;

		do{

			i2 = Select(N, rnd);

		}while(i2 == i1);

		Equal(pathnew,path,i,i1,mapping);
		Equal(pathnew,path,i+N/2,i2,mapping);

		double p = rnd->Rannyu();

		if(p < pc){

			int n = (int)rnd->Rannyu(2,l/2);
			Crossover(pathnew,i,i+N/2,l,n);

		}

	}

	//fase di mutazione genetica
	for(int i = 0; i < N; i++){

		double p = rnd->Rannyu();

		if(p < pm_pair){

			int i1 = (int)rnd->Rannyu(2,l);
			int i2 = (int)rnd->Rannyu(2,l);

			pathnew[i].Pair(i1,i2);

		}

		p = rnd->Rannyu();

		if(p < pm_shift){

			int n = (int)rnd->Rannyu(1,l/2);
			int m = (int)rnd->Rannyu(1,l-1);

			pathnew[i].Shift(n,m);

		}

		p = rnd->Rannyu();

		if(p < pm_perm){

			int m = (int)rnd->Rannyu(2,l/2);

			pathnew[i].Perm(m);


		}

		p = rnd->Rannyu();

		if(p < pm_inversion){

			int m = (int)rnd->Rannyu(2,l);

			pathnew[i].Inversion(m);

		}
	}

	//sostituisco la vecchia generazione con la nuova
	for(int i = 0; i < N; i++){

		Equal(path,pathnew,i,i,mapping);

	}

	Order(N,path); //ordino la catena

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

void Equal(Path *path1, Path *path2, int i1, int i2, map<int,City> mapping){

	int l = path2[i2].Getl();

	//copio path2 in path1
	for(int i = 0; i < l; i++){

		int c = path2[i2].GetGene(i);
		path1[i1].SetGene(i,c);

	}

	path1[i1].Setl(l);
	path1[i1].SetLenght(mapping);
}

void WriteBest(Path path, int istep){

	ofstream file;
	file.open("best.path.0",ios::app);
	double best_lenght = path.GetLenght();

	file << istep << " " << best_lenght << endl;
	file.close();

}

void WriteMean(Path *path, int istep){

	ofstream file;
	file.open("mean.path.0",ios::app);
	double sum = 0.0, sum2 = 0.0, error;

	for(int j = 0; j < 50; j++){

		double lenght = path[j].GetLenght();
		sum += lenght;
		sum2 += lenght*lenght;

	}

	sum /= (double)50;
	sum2 /= (double)50;
	error = sqrt(abs(sum2-sum*sum));
	file << istep << " " << sum << " " << error << endl;
	file.close();

}

void WriteFinal(Path path, int l, map<int,City> mapping){

	ofstream file;
	file.open("final.path.0");
	map<int,City>::iterator itr;

	for(int i = 0; i < l+1; i++){

		itr = mapping.find(path.GetGene(PBC(i,l)));
		City c = itr->second;

		double x = c.GetX();
		double y = c.GetY();

		file << i << " " << x << " " << y << endl;

	}

	file.close();

}
