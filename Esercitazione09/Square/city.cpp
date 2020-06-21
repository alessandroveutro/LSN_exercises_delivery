#include "city.h"
#include <cmath>

void City::Set(double x,double y){

	m_x = x;
	m_y = y;

}

double City::Dist(City c){

	double x = c.GetX();
	double y = c.GetY();

	double dx = x - m_x;
	double dy = y - m_y;

	return sqrt(dx*dx + dy*dy);

}
