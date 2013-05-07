/* An n-body ephemeris generation program. A Candy-Rozmus symplectic fourth order integrator is used.
To run from the command line, type: nbody [options]
options:
	switches:
		-t  :  displays the computation time. This is useful for benchmarking
		-dE :  displays the relative energy error. This is useful for accuracy checking (since dE should be 0)
		-dL :  displays the relative angular momentum error. This, too, should be 0.
		-GR :  takes general relativity into account. THIS DOESNT SEEM TO WORK!
	fields:
		time[times,t] =    : the ephemeris times. Times should be specified as a Julian Date (2455562.5 or JD2455562.5), 
						or a calender date, both in UTC.
						For a calender date, the format is year-month-day-hour-minute-second or year/month/day/hour/minute/second.
						Instead of '-' or '/', a colon could be used on the time part.
						Months could be specified by a number from 1 to 12, by its name, or by its abreviation (thus, 7, Jul, and July are all good)
						If the time part of the date is omitted, the program assumes a time of 00:00:00
						If the second isn't specified, the program assumes a time of h:m:0
						For internal purposes, the program computes using the uniform Terestial Time.
								NOTE: BUG: When the program outputs the ephemeris time, it should be the user specified time.
								NOTE: For vector mode, the program assumes input of TT, to be consistent with JPL Horizons
						If the time field is omitted, the program assumes the current time.				
		
		object[objects,obj,objs] =   :The objects for which the ephemeris is desired. Objects could be entered as names, or as the index numbers.
									The former is recommended.
									Valid objects are the eight major planets and the moon.
									If this field is omitted, the program outputs the ephemeris of all objects.
		mode = : The ephemeris output mode. The options are vector[vectors,vect] or observer[obs,ephem,ephemeris,ephemerides]
				The default mode is vector, which outputs heliocentric ecliptic (J2000) vectors.
				Observer mode outputs the equatorial RA and Dec (J2000).
		dt = : Specifies the time step in modified days (=days/k)

If everything is left blank, a help message is printed.

The initial conditions are obtained from JPL Horizons. All constants are values recommended by the International Astronomical Union.

Future work:
	allow continuous time ranges (2011 to 204000000000)
	Barycentric Coordinate Time (TCB): Investigate. Should we use that?
	Long doubles?
	
	These will facilitate longterm integration.
	
	add more integration methods
	Allow adjustible time steps
	Output orbital elements
	Fix the GR Bug
	More accurate moon ephemeris--the moon is hard to get!
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

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
double incline;

double dt=(1.0/40E0) *k;//Time step. This seems like a reasonable time step for the planets
//However, for the moon, a smaller time step (half of this) is required.

void calculate_a();   //updates the acceleration for all planets
void calculate_a_GR(); //updates the acceleration, but takes GR into account. It doesn't fully work
void dowork(double);	//integrates the planets for a specified amount of time
vector<string> split(string,char); //equivalent to the python split.
//It takes in a string, and splits it up using the character. Example: split("a,b,c",',')=[a,b,c]
int getobjnum(string); //Outputs the index of a planet/moon
string getobj(int); //Gets the name from the index
double JD(tm*); //Outputs the Julian date from a tm structure
double JD(double,double,double,double,double,double); //Outputs the Julian date from the year, month, day, hour, minute, second
double parsetime(string); //Outputs the Julian date from a string
int getmonth(string); //Returns the month number from the month name/abbreviation
void RADec(int,string&,string&,double); //Outputs the right ascension and declination for the object (specified by int) 

#include "vector.h"
//A mathematical vector class
#include "init.h"
//Initial conditions for all the planets, sun, and the moon.
#include "integrators.h"
//Defines the CRO integrator

vect totalL();  //Computes total angular momentum
double totalE(); //Computes total energy



void print_hline(){
	//Prints a line separator
	cout<<"----------------------------------------------";
	cout<<"\n";
}

void print_help(){
	cout<<"Ephemerides by David Ding\n\n";
	cout<<"The ephemerides are generated by a four-order symplectic Candy-Rozmus integration\n";
    cout<<"of the eight major planets.\n";
	cout<<"The initial position and velocity vectors are obtained from JPL Horizons\n";
	cout<<"for the date 2011 January 1, 00:00:00 UTC\n";
	cout<<"Last compiled on: "<<__DATE__;
}


void (*a)();  //a pointer to the acceleration function (GR or nonGR)

int main(int argc, char* argv[]){
    cout.precision(15); //Prints 15 digits for doubles
	vector<double> ts; //The ephemeris times
	vector<int> objs; //The objects wanted
	bool TOTIME=false; //activates timer
	bool SETTIME=false; //whether the user specfied the time or not
	bool SETPLANETS=false; //whether the user specified the queried objects
	bool GRON=false; //If general relativity is on
	bool PRINTE=false; //Whether or not to print energy
	bool PRINTL=false; //Whether or not to print angular momentum
	string mode="vectors"; //observer or vector mode
	if (argc==1){ //If no arguments are given, just help the user
		print_help();
		return 0;
	}
	
	ts.push_back(Tstart); //the starting time is in the stack of Ts, but is never used
	
	for (int i=1;i<argc;i++){ //Checks the command line arguments
		string cur_arg=string(argv[i]);
		vector<string> parts; //the lhs and rhs of an equal sign
		//Different switches
		if (cur_arg=="-t"){
			TOTIME=true; //timing
			continue;
		}
		if (cur_arg=="-GR"){
			GRON=true; //General Relativity
			continue;
		}
		if (cur_arg=="-dE"){
			PRINTE=true; //Check conservation of energy
			continue;
		}
		if (cur_arg=="-dL"){
			PRINTL=true; //Check conservation of angular momentum
			continue;
		}
		
		parts=split(cur_arg,'=');
		
		if (parts.size()!=2){//There could only be two sides of an equal sign!
			cout<<"Error in parsing commandline arguments.";
			return 1;
		}
		if (parts[0]=="time"||parts[0]=="times"||parts[0]=="t"){
			//specifying the time
			vector<string> times;
			times=split(parts[1],',');
			try{
				for (int i=0;i<times.size();i++){
					ts.push_back(parsetime(times[i]));
				}
				SETTIME=true;
			}
			catch(...){
				cout<<"Argument Error";
				return 1;
			}
		}
		else if (parts[0]=="object"||parts[0]=="objects"||parts[0]=="obj"||parts[0]=="objs"){
		//Specifying the objects
			vector<string> objects;
			objects=split(parts[1],',');
			try{
				for (int i=0;i<objects.size();i++){
					objs.push_back(getobjnum(objects[i]));
				}
				SETPLANETS=true;
			}
			catch(...){
				cout<<"Argument Error";
				return 1;
			}
		}		
		else if (parts[0]=="mode"){
		//vector or observer mode
			if (parts[1]=="vector" || parts[1]=="vect" || parts[1]=="vectors"){
				mode="vectors";
			}
			else if (parts[1]=="ephem" || parts[1]=="ephemeris" || parts[1]=="ephemerides" || parts[1]=="observer"||parts[1]=="obs") mode="observer";
			else{
				cout<<"Unrecognized mode.";
				return 1;
			}
		}
		else if (parts[0]=="dt"){
		//Sets dt
			dt=atof(parts[1].c_str());
		}
		else{
			cout<<"Unrecognized command line argument";
			return 1;
		}
	}
	//If time, planets, etc., are not set, go to default
	if (!SETTIME){
		time_t raw;
		time(&raw);
		ts.push_back(JD(gmtime(&raw)));
	}
	if (!SETPLANETS){
		for (int i=0;i<ncobjects;i++){
			objs.push_back(i);
		}
	}
	if (mode=="observer"){ //User probably specified observer time in UTC
		for (int i=1;i<ts.size();i++){
			ts[i]=TT(ts[i]);
		}
	}
	//Selects the acceleration function
	if (GRON) a=&calculate_a_GR;
	else a=&calculate_a;
	
	//We should sort ts
	
	//Initializes the planets
	initialize();
	//For timing purposes:
	clock_t start=clock();
	
	vect initL=totalL();
	double initE=totalE();
	
	cout<<"Please wait...";
	
	for (int i=1;i<ts.size();i++){
		dowork(k*ts[i]-k*ts[i-1]);
		
		cout<<"\r              ";//Erases the "Please wait..." string
		cout<<"\n";
		print_hline();
		
		cout<<"Time of ephemeris: JD"<<ts[i]<<"\n";
		for (int j=0;j<objs.size();j++){
			cout<<"\n";
			cout<<"                    "<<getobj(objs[j])<<": \n"; //Centers the object's name
			if (mode=="vectors"){
				vect rho=cobjects[objs[j]]->pos - cobjects[0]->pos;
				vect drho=cobjects[objs[j]]->v - cobjects[0]->v;
				drho=drho*k;
				cout<<"pos=\n["<<rho.x<<", "<<rho.y<<", "<<rho.z<<"] AU\n";
				cout<<"velocity=\n["<<drho.x<<", "<<drho.y<<", "<<drho.z<<"] AU/day\n";
			}
			else{
				incline=0.40909281391860255966307657332412-2.2707106390167001028956287541982e-4*(ts[0]-2451545.0)*(1.0/36525.0);
				string RAstring,Decstring;
				RADec(objs[j],RAstring,Decstring,incline);
				cout<<"RA: "<<RAstring<<"\n";
				cout<<"Dec: "<<Decstring<<"\n";
			}
		}
		cout<<"\n";
		cout<<"Please wait...";
	}
	
	cout<<"\r              "; //Erases the "Please wait..." part of the last input time
	
	if (TOTIME){ //Outputs computation time
		cout<<"\n";
		print_hline();
		cout<<"Computation Time: "<<(clock()-start)/float(CLOCKS_PER_SEC)<<" s\n";
	}
	if (PRINTE){ //Outputs relative energy error
		cout<<"\n";
		cout<<"Change in total energy: "<<(totalE()-initE)/initE<<"\n";
	}
	if (PRINTL){ //Outputs relative angular momentum error
		cout<<"\n";
		cout<<"Change in total angular momentum: ";
		print((totalL()-initL)*(1.0/initL.mag()));
	}
		
	return 0;
}

void RADec(int obj,string& RAstring, string& Decstring,double incline){
	vector<cobject> now; //Stores the current state of the solar system.
	for (int i=0;i<ncobjects;i++){
		now.push_back(*cobjects[i]);
	}
	vect rhov; //range vector
	double rho; //range
	rhov=cobjects[obj]->pos-now[3].pos;
	rho=rhov.mag();
	dowork(-rho/c); //light-time corrections;
	rhov=cobjects[obj]->pos-now[3].pos;
	/*
	rho=rhov.mag();
	for (int i=0;i<ncobjects;i++){
		*cobjects[i]=now[i];
	}
	dowork(-rho/c); //light-time corrections;
	rhov=cobjects[obj]->pos-now[3].pos;
	rho=rhov.mag();
	for (int i=0;i<ncobjects;i++){
		*cobjects[i]=now[i];
	}
	rhov=cobjects[obj]->pos-now[3].pos;
	*/
	rhov=vect(rhov.x,cos(incline)*rhov.y-rhov.z*sin(incline),rhov.y*sin(incline)+rhov.z*cos(incline)); //transforms to equatorial coordinates
	rhov=rhov*(1.00/rhov.mag()); //normalizes the range vector
	
	double dec;
	double ra;
	dec=asin(rhov.z);
	ra=acos(rhov.x/cos(dec));
	if (abs(sin(ra)*cos(dec)-rhov.y)>1E-7){ //Angle ambiguity check
		ra=2.0*pi-ra;
	}
	
	dec=180.0/pi*dec; //converts radians to degrees
	ra=12.0/pi*ra; //converts radians to hours
	
	while (ra<0) ra+=24; //RA is always positive
	
	int RAh,RAm;
	double RAs;
	RAh=int(ra);
	RAm=int((ra-(double)RAh)*60);
	RAs=(ra-(double)RAh-(double)RAm/60.0)*3600.;
	
	int Dech,Decm;
	double Decs;
	bool isneg=false; //Dec, by convention, is from -90 to 90
	if (dec<0){
		isneg=true;
		dec=-dec;
	}
	Dech=int(dec);
	Decm=int((dec-(double)Dech)*60);
	Decs=(dec-(double)Dech-(double)Decm/60.0)*3600.;
	
	char degsign=248; //degree symbol
	stringstream os1,os2;
	os1<<Dech<<degsign<<" "<<Decm<<"' "<<Decs<<"\"";
	Decstring=os1.str();
	if (isneg) Decstring="-"+Decstring;
	os2<<RAh<<" h "<<RAm<<" m "<<RAs<<" s";
	RAstring=os2.str();
	
	for (int i=0;i<ncobjects;i++){
		*cobjects[i]=now[i];
	}
}

double JD(tm* t){//Gives the Julian Date from UTC. 
//It doesn't work for dates before 1500s, since this algorithms supposes the Gregorian calander.
	double JD2000=2451544.500000; //Julian Date of Jan 1, 2000, 00:00:00
    double y=t->tm_year+1900; //tm_year is given in years since 1900
	double month=t->tm_mon+1;
    double day=t->tm_mday;
    double hour=t->tm_hour;
    double minute=t->tm_min;
	double second=t->tm_sec;
	
    double l=floor((y-2001.0)/4.0+1.0); //number of simple (Julian) leap years.
    double mlist[12]={31,28,31,30,31,30,31,31,30,31,30,31}; //Number of days in each month
	double sum=0;
	for (int i=0;i<month-1;i++){
		sum+=mlist[i];
	}
    double JDN=366.*l+365.*(y-2000.0-l)-floor((y-2001.0)/100.0)+floor((y-2001.0)/400.0)+sum+day-1.0; //The Julian Date number
    if ( (int(y)%4==0) &&((int(y)%100!=0) || (int(y)%400==0)) ){
        if (month>=3) JDN+=1; //Resolves leap years
	}
    return JD2000+JDN+hour/24.+minute/(24.*60.)+second/(24.*3600.);
}

double JD(double y, double month, double day, double hour, double minute, double second){
	double JD2000=2451544.500000;
    double l=floor((y-2001.0)/4.0+1.0);
    double mlist[12]={31,28,31,30,31,30,31,31,30,31,30,31};
	double sum=0;
	for (int i=0;i<month-1;i++){
		sum+=mlist[i];
	}
    double JDN=366.*l+365.*(y-2000.0-l)-floor((y-2001.0)/100.0)+floor((y-2001.0)/400.0)+sum+day-1.0;
    if ( (int(y)%4==0) &&((int(y)%100!=0) || (int(y)%400==0)) ){
        if (month>=3) JDN+=1;
	}
    return JD2000+JDN+hour/24.+minute/(24.*60.)+second/(24.*3600.);
}

double parsetime(string time){
	if (time.substr(0,2)=="JD") return atof(time.substr(2,time.size()-2).c_str());
		
	vector<string> parts;
	parts=split(time,'-'); //Asume users input in 2011-07-08-0-0-0
	if (parts.size()==1){ //If they didn't, doing 2011/07/08/0/0/0
		parts=split(time,'/');
	}
	if (parts.size()==1) return atof(time.c_str());
	double y=atof(parts[0].c_str());
	double month=getmonth(parts[1]); //If they enter 2011/July/7/0/0/0
	double day=atof(parts[2].c_str());
	if (parts.size()==3) return JD(y,month,day,0,0,0); //If they didn't specify time, assume midnight
	if (parts.size()==6){
		double hour=atof(parts[3].c_str());
		double minute=atof(parts[4].c_str());
		double second=atof(parts[5].c_str());
		return JD(y,month,day,hour,minute,second);
	}
	if (parts.size()!=4) throw 1; //No Idea what could go wrong :)
	vector<string> timepart;
	timepart=split(parts[3],':'); //If they enter as 2011/07/09/00:23:30 or 2011-07-09-00:23:30
	if (timepart.size()==2){ //If seconds wasn't specified, assume 0
		double hour=atof(timepart[0].c_str());
		double minute=atof(timepart[1].c_str());
		return JD(y,month,day,hour,minute,0);
	}
	if (timepart.size()==3){
		double hour=atof(timepart[0].c_str());
		double minute=atof(timepart[1].c_str());
		double second=atof(timepart[2].c_str());
		
		return JD(y,month,day,hour,minute,second);
	}
	if (timepart.size()!=1) throw 1;
	//If users input in h,m,s
	double hour=atof(split(timepart[0],'h')[0].c_str());
	double minute=atof(split(split(timepart[0],'h')[1],'m')[0].c_str());
	double second=atof(split(split(split(timepart[0],'h')[1],'m')[1],'s')[0].c_str());
	return JD(y,month,day,hour,minute,second);
	
	
}

int getmonth(string name){
	for (int i=0;i<name.size();i++){
		name[i]=(char)tolower(name[i]);
	}
	if (name=="jan" || name=="january") return 1;
	if (name=="feb" || name=="february") return 2;
	if (name=="mar" || name=="march") return 3;
	if (name=="apr" || name=="april") return 4;
	if (name=="may") return 5;
	if (name=="jun" || name=="june") return 6;
	if (name=="jul" || name=="july") return 7;
	if (name=="aug" || name=="august") return 8;
	if (name=="sep" || name=="september") return 9;
	if (name=="oct" || name=="october") return 10;
	if (name=="nov" || name=="november") return 11;
	if (name=="dec" || name=="december") return 12;
	int month;
	month=atoi(name.c_str());
	if (month<1) throw 1;
	if (month>12) throw 1;
	return month;
}

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

void calculate_a(){
	for (int j1=0;j1<ncobjects;j1++){
		cobjects[j1]->a=vect(0,0,0);
	}
	for (int j1=0; j1<ncobjects;j1++){
		for (int j2=j1+1;j2<ncobjects;j2++){
			double m1=cobjects[j1]->m;
			double m2=cobjects[j2]->m;
			vect dist=cobjects[j1]->pos-cobjects[j2]->pos;
			double magd=dist.mag();
			vect base=dist*(1.0/(magd*magd*magd));
			cobjects[j1]->a+=base*(-m2);
			cobjects[j2]->a+=base*m1;
		}
	}
}

void calculate_a_GR(){

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

void dowork(double t){
	int numtimes=int(abs(t/dt));
	dt=t/double(numtimes+1);
	numtimes=numtimes+1;
	for (int i=0;i<numtimes;i++){
		CRO_step(dt,a);
	}
}

vector<string> split(string a,char t){
	string::iterator it;
	vector<string> container;
	string cur="";
	for (it=a.begin();it!=a.end();it++){
		if (*it==t){
			container.push_back(cur); 
			cur="";
		}
		else cur+=*it;
	}
	container.push_back(cur);
	return container;
}

int getobjnum(string name){
	for (int i=0;i<name.size();i++){
		name[i]=(char)tolower(name[i]);
	}
	if (name=="sun") return 0;
	if (name=="mercury") return 1;
	if (name=="venus") return 2;
	if (name=="earth") return 3;
	if (name=="moon") return 4;
	if (name=="mars") return 5;
	if (name=="jupiter") return 6;
	if (name=="saturn") return 7;
	if (name=="uranus") return 8;
	if (name=="neptune") return 9;
	if (name=="ceres") return 10;
	if (name=="test") return 11;
	try{
		int objnum=atoi(name.c_str());
		return objnum;
	}
	catch(...){
		throw 1;
	}
}

string getobj(int o){
	if (o==0) return "Sun";
	if (o==1) return "Mercury";
	if (o==2) return "Venus";
	if (o==3) return "Earth";
	if (o==4) return "Moon";
	if (o==5) return "Mars";
	if (o==6) return "Jupiter";
	if (o==7) return "Saturn";
	if (o==8) return "Uranus";
	if (o==9) return "Neptune";
	if (o==10) return "Ceres";
	if (o==11) return "test";
	return "BAD OBJECT";
}
