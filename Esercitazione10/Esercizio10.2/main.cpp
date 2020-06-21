#include "libreria.h"
#include "mpi.h"

#define Pi 3.14152629

int main(int argc, char *argv[]){

	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();

	if(size!=4){

		if(rank == 0) cout << "Servono 4 processi, non " << size << "!!" << endl;

		return 1;

	}

	//setto il generatore di numeri casuali per ogni processo con semi differenti
	Random *rnd;
	rnd = new Random(); 
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	if (Primes.is_open()){

		for(int i = 0; i < rank; i++)

			Primes >> p1 >> p2 ;
		
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
			rnd->SetRandom(seed,p1,p2);
		}
	}

	input.close();

	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	//inizio esercizio

	int N = 100; //popolazione
	int l = 32; //numero città
	int nsteps = 100000; //numero step di simulazione

	int best[l];
	int best2[l];
	int Nmigr = 10000;

	int itag1 = 1;
	int itag2 = 2;
	int itag3 = 3;
	int itag4 = 4;

	double a = 1.0;

	map<int,City> mapping;

	//genero N città all'interno di un quadrato
	for(int i = 1; i < l+1; i++){

		int x = 0, y = 1;
		double pos[2];

		if(rank == 0){

			pos[x] = rnd->Rannyu(0,a); //coordinata x
			pos[y] = rnd->Rannyu(0,a); //coordinata y

			//double theta = rnd->Rannyu(-Pi,Pi);
			//pos[x] = a*cos(theta); //coordinata x
			//pos[y] = a*sin(theta); //coordinata y


		}

		MPI_Bcast(pos, 2, MPI_DOUBLE, 0, MPI::COMM_WORLD);

		City c;
		c.Set(pos[x],pos[y]);
		mapping.insert(pair<int,City>(i,c));

	}

	Path path[N];

	if(rank == 0){

		cout << "Simulazione algoritmo genetico" << endl;
		cout << "Numero città da visitare: " << l << endl;
		cout << "Popolazione campione: " << N << endl;
		cout << "Step evolutivi: " << nsteps << endl;
		cout << "Step migrazione: " << Nmigr << endl;
		cout << "Genero il campione iniziale di traiettorie in modo randomico..." << endl;

	}


	//genero il campione iniziale di traiettorie in modo randomico
	for(int i = 0; i < N; i++){

		path[i].SetRandom(l,rnd);
		path[i].SetLenght(mapping); //setto la lunghezza di ogni traiettoria

	}

	Order(N,path); //Ordino il mio campione per lunghezza della traiettoria

	//fase evolutiva
	for(int i = 0; i < nsteps; i++){

		WriteBest(path[0],i,rank);
		WriteMean(path,i,rank);
		Evolution(N,l,path,rnd,mapping);

		//Ogni Nmigr scambio i campioni migliori tra i processi
		if(i%Nmigr == 0 and i != 0){

			int ind[4];

			//preparo gli indici per lo scambio d'informazione tra i processi: ind[0] <-> ind[1] & ind[2] <-> ind[3]
			if(rank == 0){

				ind[0] = 0;
				ind[1] = (int)rnd->Rannyu(1,size);

				do{

					ind[2] = (int)rnd->Rannyu(1,size);

				}while(ind[2] == ind[1]);

				ind[3] = 6 - (ind[1] + ind[2]);

			}

			if(rank == 0)

				cout << "Scambio " << ind[0] << " con " << ind[1] << " e " << ind[2] << " con " << ind[3] << endl;

			//Comunico l'informazione sullo scambio a tutti i processi
			MPI_Bcast(ind, 4, MPI_INTEGER, 0, MPI::COMM_WORLD);

			//Carico l'informzaione relativa al miglior percorso di ogni processo
			for(int j = 0; j < l; j++){

				best[j] = path[0].GetGene(j);
				best2[j] = path[0].GetGene(j);

			}

			//Avvio lo scambio
			if(rank == ind[0]){

				MPI::COMM_WORLD.Isend(&best[0],l,MPI::INTEGER,ind[1],itag1);
				MPI::COMM_WORLD.Recv(&best2[0],l,MPI::INTEGER,ind[1],itag2);

			}

			if(rank == ind[1]){

				MPI::COMM_WORLD.Send(&best2[0],l,MPI::INTEGER,ind[0],itag2);
				MPI::COMM_WORLD.Recv(&best[0],l,MPI::INTEGER,ind[0],itag1);

			}

			if(rank == ind[2]){

				MPI::COMM_WORLD.Isend(&best[0],l,MPI::INTEGER,ind[3],itag3);
				MPI::COMM_WORLD.Recv(&best2[0],l,MPI::INTEGER,ind[3],itag4);

			}

			if(rank == ind[3]){

				MPI::COMM_WORLD.Send(&best2[0],l,MPI::INTEGER,ind[2],itag4);
				MPI::COMM_WORLD.Recv(&best[0],l,MPI::INTEGER,ind[2],itag3);

			}

			cout << "Processo " << rank << ": scambio terminato!" << endl;

			if(rank == ind[0] or rank == ind[2]){

				for(int j = 0; j < l; j++){

					int c = best2[j];
					path[0].SetGene(j,c);

				}
			}

			if(rank == ind[1] or rank == ind[3]){

				for(int j = 0; j < l; j++){

					int c = best[j];
					path[0].SetGene(j,c);

				}
			}
		}
	}

	if(rank == 0) cout << "Fine simulazione!" << endl;

	Order(N,path); //Ordino il mio campione per lunghezza della traiettoria
	WriteFinal(path[0],l,mapping,rank); //scrivo su file le coordinate del percorso finale
	rnd->SaveSeed(); //salvo il seme della sequenza
	MPI::Finalize();

	return 0;

}
