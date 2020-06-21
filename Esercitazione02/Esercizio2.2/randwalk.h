#ifndef __randwalk_h__
#define __randwalk_h__

#include <cmath>

#define pi 3.14152629

class RandomWalk{

	public:
		
		RandomWalk();

		void StepDiscrete(double random, double a);
		void StepContinuum(double random1, double random2, double a);
		void Reset();
		double GetR2(){return m_r2;};

	private:

		double m_x, m_y, m_z, m_r2;
		int m_step;

};



#endif
