//Initializes the planets


struct cobject{
		vect pos;
		vect pos1;
		vect v;
		vect a;
		/*cobject(double mass){
			m=mass;
		}*/
		vect fd;
		vect a1;
		float m;
};


/*
int ncobjects=10;
cobject sun(1.0);
cobject mercury(1.66019749E-7);
cobject venus(2.44806226E-6);
cobject earth(3.00374072E-6);
cobject moon(3.70050077E-8);
cobject mars(3.22742996E-7);
cobject jupiter(0.000954638698);
cobject saturn(0.000285838546);
cobject uranus(4.36664119E-5);
cobject neptune(5.15053396E-5);
*/
int ncobjects=14;
/*
cobject sun(1.0);
cobject mercury(1.0/6023600.0);
cobject venus(1.0/408523.71);
cobject earth(1.0/332946.050895);
cobject moon(3.6943037003484319109593079185365e-8);
cobject mars(1.0/3098708.0);
cobject jupiter(1.0/1047.3486);
cobject saturn(1.0/3497.898);
cobject uranus(1.0/22902.98);
cobject neptune(1.0/19412.24);
cobject ceres(4.39E-10);
*/
cobject sun,mercury,venus,earth,moon,mars,jupiter,saturn,uranus,neptune,ceres, pallas, vesta, testparticle;
cobject *cobjects[14]={&sun,&mercury,&venus,&earth,&moon,&mars,&jupiter,&saturn,&uranus,&neptune,&ceres,&pallas,&vesta,&testparticle};


void initialize(){
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
	testparticle.m=0;
	
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
	
	testparticle.pos=vect(1.0,1.0,0.0);
	testparticle.v=vect(-1.0,1.0,0.0);
	testparticle.a=vect(0,0,0);
}