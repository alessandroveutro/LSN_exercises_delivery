#include "point3D.h"

Point3D::Point3D(double x, double y, double z){

	m_x = x;
	m_y = y;
	m_z = z;

	m_r = sqrt(x*x + y*y + z*z);
	m_theta = acos(z/m_r);
	m_phi = atan(y/x);

}

void Point3D::SetCoord(double x, double y, double z){

	m_x = x;
	m_y = y;
	m_z = z;

	m_r = sqrt(x*x + y*y + z*z);
	m_theta = acos(z/m_r);
	m_phi = atan(y/x);

}

void Point3D::SetStep(double step){

	m_step = step;

}
