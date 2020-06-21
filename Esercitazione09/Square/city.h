#ifndef __city_h__
#define __city_h__

class City{

	public:
		
		void Set(double,double);
		double GetX(){return m_x;};
		double GetY(){return m_y;};
		double Dist(City c);

	private:

		double m_x, m_y;

};



#endif
