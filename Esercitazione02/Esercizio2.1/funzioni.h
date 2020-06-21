#ifndef __funzioni_h__
#define __funzioni_h__

#include <cmath>

#define pi 3.14152629

using namespace std;

class Function1D{

	public:
	
		virtual double Eval(double x) const = 0;

};

class Cos1D: public Function1D{

	public:

		double Eval(double x) const {return (pi/2)*cos((pi/2)*x);};

};

class ImportanceFun: public Function1D{

	public:

		double Eval(double x) const {return (pi/4)*cos((pi/2)*x)/(1-x);}; //funzione ottenuta invertendo la cumulativa di una distribuzione
										  //di probabilità rettilinea

};

#endif
