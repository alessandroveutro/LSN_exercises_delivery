#include "path.h"
#include "libreria.h"

void Path::SetRandom(int l, Random *rnd){

	m_genes = new int[l];
	m_l = l;

	for(int i = 0; i < l; i++){

		m_genes[i] = i+1;

	}

	int k = 0*(int)rnd->Rannyu(1,l*10);

	for(int j = 0; j < k; j++){

		int i1 = (int)rnd->Rannyu(1,l);
		int i2 = (int)rnd->Rannyu(1,l);
		Pair(i1,i2);

	}
}

void Path::SetLenght(map<int,City> mapping){

	map<int,City>::iterator itr1;
	map<int,City>::iterator itr2;
	double d = 0.0;
	int l = m_l;

	for(int i = 0; i < l; i++){

		itr1 = mapping.find(m_genes[PBC(i,l)]);
		itr2 = mapping.find(m_genes[PBC(i+1,l)]);

		City c1 = itr1->second;
		City c2 = itr2->second;

		d += c1.Dist(c2);

	}

	m_lenght = d;
}

void Path::SetGene(int n, int c){

	m_genes[n] = c;

}

void Path::Setl(int l){

	m_l = l;

}

void Path::Pair(int i1, int i2){

	double k = m_genes[i1];
	m_genes[i1] = m_genes[i2];
	m_genes[i2] = k;

}

void Path::Shift(int n, int m){

	//1 < n < l/2
	//1 < m < l-1

	int genes[m];

	for(int i = 1; i < m+1; i++){

		int k = m_genes[i];
		genes[i] = k;
		
		if(k + n > m_l)

			m_genes[i] = k + n + 1 - m_l;

		else

			m_genes[i] = k + n;

	}

	int ind = m+1;

	do{

		int count = 0;

		for(int i = 1; i < ind; i++)

			if(m_genes[ind] == m_genes[i])

				count++;

		if(count == 0){

			ind++;

		}else if(count == 1){

			for(int j = 1; j < m_l; j++){

				count = 0;
	
				for(int i = 1; i < ind; i++)

					if(genes[j] == m_genes[i])

						count++;

				if(count == 0){

					m_genes[ind] = genes[j];
					ind++;
					break;

				}
			}
		}

	}while(ind < m_l);
}

void Path::Perm(int m){

	//2 < m < l/2

	for(int i = 1; i < m; i++){

		int k = m_genes[i];
		m_genes[i] = m_genes[m_l-i];
		m_genes[m_l-i] = k;

	}

}

void Path::Inversion(int m){

	//2 < m < l

	for(int i = 1; i < m; i++){

		int k = m_genes[i];
		m_genes[i] = m_genes[m-i];
		m_genes[m-i] = k;

	}

}
