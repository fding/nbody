import nbody,time,math



def sqerror(v1,v2):
    return math.sqrt(sum([(v1[i]-v2[i])**2 for i in range(3)]))*149597871


start=time.clock()
sys=nbody.SolarSystem(2459214.9)
end=time.clock()
#_titan_triton_plutoBC
f=open('errors_planets_and_galileanmoons_saturnian_triton_uranianmoons_asteroids.txt','w')
f.write("Time: 2459214.9")
f.write("\nSimulation time: "+str(end-start)+" seconds")

pos=[0]*len(sys)
pos[0]=(-6.647122354508441E-03, 6.007351268184804E-03, 1.055511212196090E-04)
pos[1]=(2.195061271263852E-01,-3.618867074836926E-01,-5.070256315775003E-02)
pos[2]=(-4.624687396981217E-01,-5.562863018697876E-01,1.869270037060782E-02)
pos[3]=(-1.754388146172777E-01,9.746833297341465E-01,5.887240751783555E-05)
pos[4]=(-1.765198314300658E-01,9.770405480026660E-01,1.954461691774774E-04)
pos[5]=(6.215534430239492E-01,1.377492273505372E+00,1.343595416115652E-02)
pos[6]=(3.031250844135190E+00,-4.085073313290786E+00,-5.087071789243780E-02)
pos[7]=(5.481098580161068E+00,-8.337949724855315E+00,-7.323611545819925E-02)
pos[8]=(1.534487597728014E+01,1.246589390197973E+01,-1.524962389740150E-01)
pos[9]=(2.945358676754011E+01,-5.227714696709961E+00,-5.711243647369795E-01)
for i in range(10):
    f.write("\nError for "+sys[i].name+": "+str(sqerror(pos[i],sys[i].position))+' km')

f.flush()
f.close()
