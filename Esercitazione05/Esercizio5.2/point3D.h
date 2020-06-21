#ifndef __point3D_h__
#define __point3D_h__

#include <cmath>

#define pi 3.14152629

class Point3D{

	public:
		
		Point3D(double x, double y, double z);

		void SetCoord(double x, double y, double z);
		void SetStep(double step);

		double GetX(){return m_x;};
		double GetY(){return m_y;};
		double GetZ(){return m_z;};
		double GetR(){return m_r;};
		double GetTheta(){return m_theta;};
		double GetPhi(){return m_phi;};
		double GetStep(){return m_step;};

	private:

		double m_x, m_y, m_z;
		double m_r, m_theta, m_phi;
		double m_step;

};



#endif
