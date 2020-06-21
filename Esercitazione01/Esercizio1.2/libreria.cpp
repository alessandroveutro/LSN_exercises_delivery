#include "libreria.h"

void StandardDice(double *r, Random rnd, int dim){

	for(int i = 0; i < dim; i++)

		r[i] = (int)rnd.Rannyu(1,7); //genero un numero (pseudo)casuale intero tra 1 e 6, i.e. le facce di un dado

}

void ExpDice(double *r, Random rnd, int dim, double lambda){

	for(int i = 0; i < dim; i++)

		r[i] = rnd.Exp(lambda); 

}

void CauchyDice(double *r, Random rnd, int dim, double mu, double gamma){

	for(int i = 0; i < dim; i++)

		r[i] = rnd.Cauchy(mu, gamma); 

}
