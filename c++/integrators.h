//integrators

void CRO_step(register double mydt,void (*a)()){
	long double macr_a[4] = {0.5153528374311229364, -0.085782019412973646,0.4415830236164665242, 0.1288461583653841854};
	long double macr_b[4] = {0.1344961992774310892, -0.2248198030794208058, 0.7563200005156682911, 0.3340036032863214255};
	for (int i=0;i<4;i++){
            a();
            for (int j=0;j<ncobjects;j++){
                cobjects[j]->v +=  cobjects[j]->a * mydt*macr_b[i];
                cobjects[j]->pos += cobjects[j]->v  * mydt*macr_a[i]; 
			}
	} //We should really expand the loop for efficiency
}



/*
void verlet(double t){
	initialize();
	Checks for momentum conservation
    vect p;
    print(p);

	if (k*t<Tstart) dt=-dt;
	if (k*t==Tstart) return;
	int eulern=1;
	double idt=dt/eulern;
	for (int i=0;i<eulern;i++){
		for (int j1=0; j1<ncobjects;j1++){
			for (int j2=j1+1;j2<ncobjects;j2++){
				vect dist=cobjects[j1]->pos-cobjects[j2]->pos;
				double magd=dist.mag();
				vect base=dist*(1.0/(magd*magd*magd));
				double h=cross(dist,cobjects[j1]->v-cobjects[j2]->v).mag();
			    vect base1=dist*(3.0*h*h*(cobjects[j1]->m+cobjects[j2]->m)/(c2*magd*magd*magd*magd*magd));
				cobjects[j1]->a+=(base*(-cobjects[j2]->m)-base1*(cobjects[j2]->m));
				cobjects[j2]->a+=(base*cobjects[j1]->m+base1*cobjects[j1]->m);

				cobjects[j1]->a+=base*(-cobjects[j2]->m);
				cobjects[j2]->a+=base*cobjects[j1]->m;
			}
			cobjects[j1]->v+=cobjects[j1]->a*idt;
			cobjects[j1]->pos+=cobjects[j1]->v*idt-cobjects[j1]->a*idt*(0.5*idt);
			cobjects[j1]->a=vect(0,0,0);
		}
	}
	t_interval+=dt;
	double dt2=dt*dt;
	
	while (fabs(t_interval)<fabs(runt_interval)){
		t_interval+=dt;
		for (int j1=0; j1<ncobjects;j1++){
			for (int j2=j1+1;j2<ncobjects;j2++){
                double m1=cobjects[j1]->m;
                double m2=cobjects[j2]->m;
				vect dist=cobjects[j1]->pos-cobjects[j2]->pos;
				double magd=dist.mag();
				vect base=dist*(1.0/(magd*magd*magd));
				cobjects[j1]->a+=base*(-cobjects[j2]->m);
				cobjects[j2]->a+=base*cobjects[j1]->m;
				double h=cross(dist,cobjects[j1]->v-cobjects[j2]->v).mag();
                vect base1=dist*(3.0*h*h*(m1+m2)/(c2*pow(magd,5)));
				cobjects[j1]->a+=(base*(-m2)-base1*m2);
				cobjects[j2]->a+=(base*m1+base1*m1);
				
				
				Fifth order correction  
				double rrdot=dot(dist,cobjects[j1]->v-cobjects[j2]->v);
				cobjects[j1]->fd+=dist*((-2*m2*m2)/pow(magd,6)+3*m2*rrdot/pow(magd,5)-15*m2*pow(rrdot,2)/pow(magd,7));
				cobjects[j1]->fd+=(cobjects[j1]->v-cobjects[j2]->v)*6.0*m2*(rrdot/pow(magd,5));
				cobjects[j2]->fd+=dist*((2*m1*m1)/pow(magd,6)-3*m1*rrdot/pow(magd,5)+15*m1*pow(rrdot,2)/pow(magd,7));
				cobjects[j2]->fd+=(cobjects[j2]->v-cobjects[j1]->v)*6.0*m1*(rrdot/pow(magd,5));
			}
			vect temp=cobjects[j1]->pos*2.0-cobjects[j1]->pos1+cobjects[j1]->a*dt2;//+cobjects[j1]->fd*(0.08333333*dt2*dt2);
			
            cobjects[j1]->v=(temp-cobjects[j1]->pos1)*(0.5/dt);
			cobjects[j1]->pos1=cobjects[j1]->pos;
			cobjects[j1]->pos=temp;
			cobjects[j1]->a=vect(0,0,0);
		}
	}
    Checks for momentum conservation
	p=vect(0,0,0);
	for (int i=0;i<ncobjects;i++){
        p+=cobjects[i]->v*cobjects[i]->m;
    }
    print(p); 
}*/
