'''
Todo:
1. Friendlier time system, besides just accepting Julian Dates. For instance, maybe create a Time class?
2. Vector class?
3. SolarSystem function should take more sets of initial conditions to save time.
'''

import sys
import datetime as pyt

if sys.platform=='win32':
    import windows.nbody as cnbody
elif sys.platform=='linux2' or sys.platform=='linux':
    print 'Linux is not supported yet, sorry'
    sys.exit()
elif sys.platform=='cygwin':
    print 'Cygwin is not supported yet, sorry'
    sys.exit()
elif sys.platform=='darwin':
    print 'Mac OS X is not supported yet, sorry'
    sys.exit()
elif sys.platform=='os2':
    print 'OS/2 is not supported yet, sorry'
    sys.exit()
elif sys.platform=='os2emx':
    print 'OS/2 EMX is not supported yet, sorry'
    sys.exit()
elif sys.platform=='riscos':
    print 'RiscOS is not supported yet, sorry'
    sys.exit()
elif sys.platform=='atheos':
    print 'AtheOS is not supported yet, sorry'
    sys.exit()

_K = 0.01720209895

def _vect_to_tuple(vect):
    return (vect.x,vect.y,vect.z)

class CObject:
    def __init__(self,mass=0.0, position=(0,0,0), velocity=(0,0,0),name="<celestial body>"):
        self._cobject = cnbody.cobject()
        self.mass=mass
        self.position=position
        self.velocity=velocity
        self.name=name
    def __setattr__(self, name, value):
        self.__dict__[name]=value
        if name=="mass":
            self._cobject.m=value
        if name=="position":
            self._cobject.pos=cnbody.vect(value[0],value[1],value[2])
        if name=="velocity":
            self._cobject.v=cnbody.vect(value[0]/_K,value[1]/_K,value[2]/_K)
    def _update(self):
        self.mass=self._cobject.m
        self.velocity=_vect_to_tuple(self._cobject.v*_K)
        self.position= _vect_to_tuple(self._cobject.pos)
    def __str__(self):
        return "[Object "+self.name+": Mass="+str(self.mass)+" Position="+str(self.position)+" Velocity="+str(self.velocity)+"]"

class System:
    def __init__(self,objects=[], time=pyt.date.today()):
        self._System= cnbody.System()
        self.objects=objects
        self.time = time #Re-examine this line, to allow for more general time inputs
        self._System.settime(self.time)
        for obj in objects:
            self._System.addobject(obj._cobject)
    def simulate(self,time=pyt.date.today(), dt=0.05):
        self._System.setdt(dt)
        self._System.runtotime(time) #Re-examine this line, to allow for more general time inputs
        for obj in self.objects:
            obj._update()
        self.time=self._System.gettime()
    def __str__(self):
        returnstring= "System Time="+str(self.time)
        for obj in self.objects:
            returnstring+="\n"+str(obj)
        return returnstring

def SolarSystem(time=2455562.500000000):
    execfile("init_conditions/2455562.5.py",globals())
    sys=System(time=2455562.500000000,objects=[Sun,Mercury,Venus,Earth,Moon,Mars,Jupiter,Saturn,Uranus,Neptune,Ceres,Pallas,Vesta])
    sys.simulate(time=time)
    return sys
