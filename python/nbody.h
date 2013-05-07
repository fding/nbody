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



cobject sun,mercury,venus,earth,moon,mars,jupiter,saturn,uranus,neptune,ceres, pallas, vesta;
	
System* SolarSystem(){
	cobject* ssobjects[13];
	sun.m=1.0;
	mercury.m=1.0/6023600.0;
	venus.m=1.0/408523.71;
	earth.m=1.0/332946.050895;
	moon.m=3.6943037003484319109593079185365e-8;
	mars.m=1.0/3098708.0;
	jupiter.m=1.0/1047.3486;
	saturn.m=1.0/3497.898;
	uranus.m=1.0/22902.98;
	neptune.m=1.0/19412.24;
	ceres.m=4.7E-10;
	pallas.m=1.0E-10;
	vesta.m=1.3E-10;

	sun.pos= vect(-4.173780072034275E-03,7.099593014214438E-04,1.604504857505994E-05);
	sun.v=vect(8.230492550577723E-07,-6.168298835870542E-06,-8.099453361184341E-09)*(1.0/k);
	sun.a=vect(0,0,0);

	mercury.pos=vect(-3.254274813047346E-01,1.470377149611438E-01,4.144867298810978E-02);
	mercury.v=vect(-1.745330858881957E-02,-2.441610498768040E-02,-3.928875849872461E-04)*(1.0/k);
	mercury.a=vect(0,0,0);

	venus.pos=vect(-5.422830038348582E-01,4.753132014458776E-01,3.757371712506968E-02);
	venus.v=vect(-1.345295870164601E-02,-1.527996051896081E-02,5.672054215939561E-04)*(1.0/k);
	venus.a=vect(0,0,0);

	earth.pos=vect(-1.757691570699245E-01,9.689784107710354E-01,-8.071357286641453E-06);
	earth.v=vect(-1.722543139856719E-02,-3.069797666532440E-03,-4.254847485630660E-07)*(1.0/k);
	earth.a=vect(0,0,0);

	moon.pos=vect(-1.770707371941026E-01,9.668060600096702E-01,-1.373909314834243E-04);
	moon.v=vect(-1.672410611695570E-02,-3.394944511095675E-03,4.453219362728522E-05)*(1.0/k);
	moon.a=vect(0,0,0);

	mars.pos=vect(5.680571685512422E-01,-1.290276310379983E+00,-4.108702446348692E-02);
	mars.v=vect(1.332175320659049E-02,6.865622414186563E-03,-1.831125709754552E-04)*(1.0/k);
	mars.a=vect(0,0,0);

	jupiter.pos=vect(4.901953649524238E+00,6.492361425410386E-01,-1.124667705413734E-01);
	jupiter.v=vect(-1.081844011900388E-03,7.839399858800254E-03,-8.404735693140495E-06)*(1.0/k);
	jupiter.a=vect(0,0,0);

	saturn.pos=vect(-9.415804158693506E+00,-1.768796125725446E+00,4.054587777717821E-01);
	saturn.v=vect(7.310506583016966E-04,-5.494644595066728E-03,6.629310885757098E-05)*(1.0/k);
	saturn.a=vect(0,0,0);

	uranus.pos=vect(2.008289298043070E+01, -1.626990475698014E-01,-2.607906257647571E-01);
	uranus.v=vect(2.987279739340398E-06,3.749614174715089E-03,1.393893130105414E-05)*(1.0/k);
	uranus.a=vect(0,0,0);

	neptune.pos=vect(2.543242295413326E+01,-1.592871746992481E+01,-2.580858664951978E-01);
	neptune.v=vect(1.645904065912284E-03,2.679121207007654E-03,-9.328009914558260E-05)*(1.0/k);
	neptune.a=vect(0,0,0);

	ceres.pos=vect(1.678042856488875E+00,-2.398633539551728E+00,-3.848350054846178E-01);
	ceres.v=vect(7.953867896114812E-03,5.313406707674584E-03 ,-1.299702261339850E-03)*(1.0/k);
	ceres.a=vect(0,0,0);

	pallas.pos=vect(1.021955602810206E-01,-2.682258151022251E+00,1.845303594657398E+00);
	pallas.v=vect(8.571470711239143E-03,-1.196471935945572E-03,1.088113157993274E-04);
	pallas.a=vect(0,0,0);

	vesta.pos=vect(-9.417630876151264E-01,-1.932215041162336E+00,1.720699960471735E-01);
	vesta.v=vect(1.097668925917006E-02,-5.244970192356581E-03,-1.175969509428757E-03);
	vesta.a=vect(0,0,0);
	
	ssobjects[0]=&sun;
	ssobjects[1]=&mercury;
	ssobjects[2]=&venus;
	ssobjects[3]=&earth;
	ssobjects[4]=&moon;
	ssobjects[5]=&mars;
	ssobjects[6]=&jupiter;
	ssobjects[7]=&saturn;
	ssobjects[8]=&uranus;
	ssobjects[9]=&neptune;
	ssobjects[10]=&ceres;
	ssobjects[11]=&pallas;
	ssobjects[12]=&vesta;
	System *solarSystem=new System(vector<cobject*>(ssobjects,ssobjects+13));
	solarSystem->settime(Tstart);
	return solarSystem;
}



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
	if (numtimes>100000) numtimes=100000;
	dt=(k*t-k*this->time)/double(numtimes+1);
	numtimes=numtimes+1;
	//if (GRON) a=&System::calculate_a_GR;
	//else a=&System::calculate_a;
	for (int i=0;i<numtimes;i++){
		CRO_step(dt);
	}
	this->time=t;
}


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