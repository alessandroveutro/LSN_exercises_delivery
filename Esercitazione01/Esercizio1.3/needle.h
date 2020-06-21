#ifndef __needle_h__
#define __needle_h__

class Needle{

	public:
		
		Needle();

		void Throw(double yc, double length, double theta);
		void CheckHit(double distance);
		double GetHit(){return m_hit;};

	private:

		double m_yc, m_ymin, m_ymax;
		int m_hit;

};



#endif
