#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <sstream>

#ifndef NBODY
#define NBODY
using namespace std;

//Initializes the planets

#include <vector>

#include "vector.h"


// Not sure if needed
double TT(double t){
	//Takes in a UTC Julian date and outputs a Terestial Date Julian Date
	//Does not take leap seconds into account
    double T=(t- 2451545.0)/36525. ;
    double ET_UT=64.184+59.0* T-51.2*T*T - 67.1 * T*T*T - 16.4* T*T*T*T;
    return t+ET_UT/86400.0;
}


const double k=0.01720209895;  //Gaussian gravitational constant
const double Tstart=TT(2455562.500000000); //2011-January 1, 0:0:0 UTC
const double c=173.144483/k; //Speed of light in AU/modified day
const double c2=c*c; //Square of the speed of light, used in GR calculations
const double pi=3.1415926535897932384626433832795;

//double dt=(1.0/40E0) *k;//Time step. This seems like a reasonable time step for the planets
//However, for the moon, a smaller time step (half of this) is required.


class cobject{
	public:
		vect* getpos();
		vect pos;
		vect pos1;
		vect v;
		vect a;
		vect fd;
		vect a1;
		double m;
};

vect* cobject::getpos(){
	return &pos;
}

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

void System::addobject(cobject* newobject){
	cobjects.push_back(newobject);
}
cobject* System::getobject(int i){
	return cobjects[i];
}
void System::settime(double t){
	this->time=t;
}
double System::gettime(){
	return this->time;
}
void System::setdt(double dt){
	this->dt=dt;
}
System::System(vector<cobject*> list){
	this->cobjects=list;
	this->dt=(1.0/20E0) *k;
}
System::System(){
	this->dt=(1.0/20E0) *k;
}

// Candy-Rozmus integrator: 4-th order symplectic
void System::CRO_step(register double mydt){
	long double macr_a[4] = {0.5153528374311229364, -0.085782019412973646,0.4415830236164665242, 0.1288461583653841854};
	long double macr_b[4] = {0.1344961992774310892, -0.2248198030794208058, 0.7563200005156682911, 0.3340036032863214255};
	for (int i=0;i<4;i++){
            this->calculate_a();
			vector<cobject*>::iterator it;
            for (it=this->cobjects.begin();it!= this->cobjects.end();it++){
                (*it)->v +=  (*it)->a * mydt*macr_b[i];
                (*it)->pos += (*it)->v  * mydt*macr_a[i]; 
			}
	} //We should really expand the loop for efficiency
}


void System::runtotime(double t){
	int numtimes=int(abs(k*t-k*this->time)/this->dt);
	// Ensures that dt goes into t an integer amount of times
	// Also, we must convert to modified time for the simulation
	dt=(k*t-k*this->time)/double(numtimes+1);
	numtimes=numtimes+1;
	//if (GRON) a=&System::calculate_a_GR;
	//else a=&System::calculate_a;
	for (int i=0;i<numtimes;i++){
		CRO_step(dt);
	}
	this->time=t;
}

// Calculates the acceleration due to Newtonian gravity
void System::calculate_a(){
	vector<cobject*>::iterator it, it1;
	for (it=this->cobjects.begin();it!=this->cobjects.end();it++){
		(*it)->a=vect(0,0,0);
	}
	for (it1=this->cobjects.begin();it1!=this->cobjects.end();it1++){
		vector<cobject*>::iterator it2;
		for (it2=it1+1;it2!=this->cobjects.end();it2++){
			double m1=(*it1)->m;
			double m2=(*it2)->m;
			vect dist=(*it1)->pos-(*it2)->pos;
			double magd=dist.mag();
			vect base=dist*(1.0/(magd*magd*magd));
			(*it1)->a+=base*(-m2);
			(*it2)->a+=base*m1;
		}
	}
}

/*
void System::calculate_a_GR(){

	const double gamma=1; //For other theories of relativity
	const double beta=1;//For other theories of gravity
	calculate_a();
	double table1[ncobjects];
	double disttable[ncobjects][ncobjects];
	for (int j1=0; j1<ncobjects; j1++){
		table1[j1]=0;
		cobjects[j1]->a1=vect(0,0,0);
		for (int j2=0; j2<ncobjects; j2++)
		{
			if (j2==j1) continue;
			vect dist=cobjects[j1]->pos-cobjects[j2]->pos;
			double magd=dist.mag();
			disttable[j1][j2]=magd;
			table1[j1]+=cobjects[j2]->m/magd;
		}
	}
	for (int j1=0; j1<ncobjects;j1++){
		for (int j2=0;j2<ncobjects;j2++){
			if (j1 == j2) continue;
			double m1=cobjects[j1]->m;
			double m2=cobjects[j2]->m;
			vect dist=cobjects[j1]->pos-cobjects[j2]->pos;
			double v1=cobjects[j1]->v.mag();
			double v2=cobjects[j2]->v.mag();
			double magd=disttable[j1][j2];
			double imagd3=1.0/(magd*magd*magd);
			double vdot=dot(cobjects[j1]->v,cobjects[j2]->v);
			double dot1=dot(dist,cobjects[j2]->v);
			double dot2=dot(dist,cobjects[j2]->a);
			double base=imagd3*cobjects[j2]->m;
			double r_series1,r_series2; //Relativity;
			r_series1=1.0-2.0*(beta+gamma)/c2*table1[j1]-(2.0*beta-1)/c2*table1[j2]+gamma*v1*v1/c2+(1.0+gamma)*v2*v2/c2-2.0*(1.0+gamma)/c2*vdot-1.5/c2*(dot1*dot1/magd/magd)-0.5/c2*dot2;
			r_series2=dot(dist,(cobjects[j1]->v-cobjects[j2]->v)*(1+2*gamma)+cobjects[j1]->v);
			cobjects[j1]->a1+=dist*(-r_series1*base)+(cobjects[j1]->v-cobjects[j2]->v)*(r_series2*1.0/c2*base);
			cobjects[j1]->a1+=cobjects[j2]->a*((1.5+2.0*gamma)/c2*cobjects[j2]->m/magd);
		}
		cobjects[j1]->a=cobjects[j1]->a1;
	}
}

*/
/*

double totalE(){
	double E=0;
	for (int i=0;i<ncobjects;i++){
		for (int j=i+1;j<ncobjects;j++){
			E+= -(cobjects[i]->m*cobjects[j]->m)/(cobjects[i]->pos-cobjects[j]->pos).mag();
		}
		E+=0.5*cobjects[i]->m*pow((cobjects[i]->v).mag(),2);
	}
	return E;
}

vect totalL(){
	vect L;
	for (int i=0;i<ncobjects;i++){
		L+=cross(cobjects[i]->pos,cobjects[i]->v)*cobjects[i]->m;
	}
	return L;
}

double TT(double t){
	//Takes in a UTC Julian date and outputs a Terestial Date Julian Date
	//Does not take leap seconds into account
    double T=(t- 2451545.0)/36525. ;
    double ET_UT=64.184+59.0* T-51.2*T*T - 67.1 * T*T*T - 16.4* T*T*T*T;
    return t+ET_UT/86400.0;
}

*/

#endif