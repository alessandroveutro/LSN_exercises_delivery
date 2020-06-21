#include "needle.h"
#include <cmath>

Needle::Needle(){

	m_yc = 0;
	m_ymin = 0;
	m_ymax = 0;
	m_hit = 0;

}

void Needle::Throw(double yc, double length, double theta){

	m_yc = yc;

	if(theta > 0){

		m_ymin = m_yc - (length/2)*sin(theta);
		m_ymax = m_yc + (length/2)*sin(theta);

	}else{

		m_ymin = m_yc + (length/2)*sin(theta);
		m_ymax = m_yc - (length/2)*sin(theta);

	}
	
	m_hit = 0;

}

void Needle::CheckHit(double distance){

	for(int j = 0; m_ymax > (double)j*distance; j++)

		if(m_ymin < (double)distance*j)

			m_hit = 1;

	if(m_hit != 1)

		m_hit = 0;

}
