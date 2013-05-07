nbody
=====

Fast n-body gravitational simulation code, written in  C++ with ports to other slower languages.

Currently supported languages: C++, Python.

To use this code for Python, build the C++ files as:

swig -c++ -python nbody.i
gpp nbody_wrap.cxx.

Alternatively, if you are running Windows, just copy the precompiled python module in the windows folder 
to your working directory.

Why use this code?
=================

1) N-body simulations written entirely in Python are way too slow. They can be used to generate pretty models,
but to generate accurate data (arcsecond-level accuracy years into the future), the computation time needed is too big.

2) N-body simulations written in C++ is fast, but C++ isn't known to be an interactive language. 
If you don't need interactivity, just use the C++ code (or even just compile the C++ ephemeris program). On the other hand,
if you feel like you want to play around with the simulations, integrate it with your other code, save days of debugging time,
avoid re-compiling your code every few-minutes, etc., C++ is not a great language.

This Python (and perhaps other languages in the future) port allows you to interact with the simulations while performing
the brute work in fast compiled and optimized C++. 

Methodology
===========

This program uses Newtonian mechanics and a four-order symplectic Candy-Rozmus integration
(a symplectic algorithm guarantees exact conservation of energy and angular momentum).
The initial conditions are obtained from JPL Horizons, ahd constants (like masses, gravitational constant) are those
recommended by the International Astronomical Union. The program currently does not take into account effects like general
relativity, the non-spherical shapes of celestial objects, tidal effects on Earth, etc. It also does not take the
500 asteroids used by JPL Horizons into accound in its model of the Solar System.


Usage
=====

|Ephemeris program|
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


|Python|

You can import the module by
  import nbody
The module contains a cobject class, a System class, and a vector (vect) class.
cobject is used to store information about a celestial object -- its mass, position, and velocity, stored as
cobject.m, cobject.pos, and cobject.v.
The vect class is used to store vectors. You can access the coordinates of the vector via vect.x, vect.y,
and vect.z. You can add, subtract, multiply by scalars, etc., in the obvious way.
The System class is used to store a system of celestial objects. You can add an object to the system by calling
System.addobject(cobject). You can access the i-th member of the System by calling System.getobject(i). 
You can set the time-state of the system by calling System.settime(), and get the time by calling System.gettime().
Finally, you can call System.runtotime(t) to iterate the system to time t.

An additional useful function is SolarSystem(), which returns the solar system on January 1, 2011. The solar system
includes, in order,
Sun, Mercury, Venus, Earth, MOON, Mars, Jupiter, Saturn, Uranus, Neptune, Ceres, Pallas, Vesta.
(Pluto will be included as a feature sometime later).
Note that Mars is not SolarSystem().getobject(4); it is, rather, SolarSystem().getobject(4).

All times must be specified in Julian dates. Coordinates should be using the ecliptic coordinate system. Distances
are measured in AU, and masses in solar masses. The origin is the solar system barycenter.

Example
=======

This code will give the x, y, z location of Venus on January 1, 2014:

import nbody
v=nbody.SolarSystem().runtotime(2456659).getobject(2).pos
print (v.x,v.y,v.z)

This code will add a test particle at location (1,1) and velocity (-1,1) on January 1, 2011, and find the location
of that particle in January 1, 2014:

import nbody
solarsystem=nbody.SolarSystem()
testparticle=cobject()
testparticle.m=0
testparticle.pos=vect(1,1,0)
testparticle.v=vect(-1,1,0)
solarsystem.addobject(testparticle)
v=nbody.SolarSystem().runtotime(2456659).getobject(13).pos
print (v.x,v.y,v.z)

Bugs, deficiencies, and Future Plans
===================================

Yeah, this version of the code probably has a few bugs, and the interface with Python is rather ugly.
I plan to:
1. Create a nicer and more intuitive wrapper, maybe along the lines of:

import nbody
SolarSystem.simulate(t=today())
print Venus.position

and

import nbody
testparticle=cobject(mass=0, position=(1,1,0), velocity=(-1,1,0))
mysystem = SolarSystem.add(testparticle,time=today())
mysystem.simulate(t, dt=0.1)
print testparticle.position

2. Compile binaries for Linux and Macs.
3. Fix any bugs that may come up.



