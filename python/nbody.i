%module nbody


%{
#include "vector.h"
#include "nbody.h"
%}

class System{
	public:
		void addobject(cobject* newobject);
		cobject* getobject(int i); //returns the i-th object
		void runtotime(double); //runs simulation to time t
		void settime(double);
		double gettime();
		void setdt(double);
		void calculate_a();
		System(vector<cobject*> list);
		System();
	private:
		vector<cobject*> cobjects;
		void CRO_step(register double mydt);
		double time;
		double dt;
};


class vect{
	public:
		double x;
		double y;
		double z;
		vect();
		vect(double xin, double yin, double zin);
		bool operator==(const vect &rhs);
		vect & operator+=(const vect &rhs);
		vect & operator-=(const vect &rhs);
		vect & operator*=(const double &rhs);
		vect operator+(const vect &other);
		vect operator-(const vect &other);
		vect operator*(const double &other);
		bool operator!=(const vect& other);
		double mag();
};

class cobject{
	public:
		vect* getpos();
		vect pos;
		vect pos1;
		vect v;
		vect a;
		vect fd;
		vect a1;
		float m;
};


//extern System* SolarSystem();

extern double dot(const vect& a, const vect& b);

extern inline vect cross(const vect& a, const vect& b);
