#include "randwalk.h"

RandomWalk::RandomWalk(){

	m_x = 0;
	m_y = 0;
	m_z = 0;
	m_r2 = 0;
	m_step = 0;

}

void RandomWalk::StepDiscrete(double random, double a){

	//random deve essere un numero intero distribuito uniformemente tra 0 e 5 (6 numeri che indicano le 6 possibili direzioni)

	if(random == 0){

		m_x += a;

	}else if(random == 1){

		m_x -= a;

	}else if(random == 2){

		m_y += a;

	}else if(random == 3){

		m_y -= a;

	}else if(random == 4){

		m_z += a;

	}else if(random == 5){

		m_z -= a;

	}

	m_r2 = m_x*m_x + m_y*m_y + m_z*m_z;
	m_step ++;

}

void RandomWalk::StepContinuum(double random1, double random2, double a){

	//random1 e random 2 devono essere due numeri estratti uniformemente, rispettivamente dagli intervalli (0,2*pi) e (0,1)
	//permettono di estrapolare la posizione angolare di un punto distribuito uniformemente all'interno dell'angolo solido
	double phi = random1;
	double s = random2;
	double theta = acos(1-2*s);

	m_x += a*sin(theta)*cos(phi);
	m_y += a*sin(theta)*sin(phi);
	m_z += a*cos(theta);
	m_r2 = m_x*m_x + m_y*m_y + m_z*m_z;
	m_step ++;

}

void RandomWalk::Reset(){

	//resetto il random walk
	m_x = 0;
	m_y = 0;
	m_z = 0;
	m_r2 = 0;
	m_step = 0;

}


