#ifndef __path_h__
#define __path_h__

#include "city.h"
#include "random.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <map>
#include <iterator>

using namespace std;

class Path{

	public:
		
		void SetRandom(int l, Random *rnd);
		void SetLenght(map<int,City> mapping);
		void SetGene(int n, int c);
		void Setl(int l);
		double GetGene(int n){return m_genes[n];};
		double GetLenght(){return m_lenght;};
		int Getl(){return m_l;};

		void Pair(int i1, int i2);
		void Shift(int n, int m);
		void Perm(int m);
		void Inversion(int m);

	private:

		int *m_genes;
		double m_lenght;
		int m_l;

};

#endif
