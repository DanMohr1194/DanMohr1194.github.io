#doing calculations based on Kepler's laws with velocities, positions, and masses of the bodies

import math
import time
import sys


#Global variables
grav_const=6.674e-11 #newton's gravitational constant (N(m/kg)^2)

#the body in orbit
class Satellite:
    #parameters for the satellite, orbits host
    def __init__(self, rad_init, vel_init, mass):
        self._rad_init=rad_init #initial radius from host in +x direction, meters
        self._mass=mass #mass, kg
        self._vel_init=vel_init #initial velocity, m/s, always vertical (on the screen)

#the body that's being orbited around
class Host:
    #host body, orbited by satellite
    def __init__(self, mass):
        self._mass=mass #mass, kg
    
    #escape velocity from host, m/s
    def escape_vel(self, rad):
        return math.sqrt(2 * grav_const * self._host._mass / rad)

#the system including the host and satellite
class Orbital_System:
    def __init__(self, sat, host):
        self._sat=sat
        self._host=host

        #unchanging orbital elements
        self._grav_pams = grav_const*(self._sat._mass+self._host._mass) #sum of standard gravitational parameters
        self._spec_orb_energ = (self._sat._vel_init ** 2) / (2) - self._grav_pams / sat._rad_init #specific orbital energy, J/kg
        self._ang_mom = self._sat._rad_init * self._sat._vel_init  # specific relative angular momentum, m^2/s
        self._eccen = math.sqrt((1 + (((2 * self._spec_orb_energ) * self._ang_mom ** 2) / (self._grav_pams ** 2)))) #eccentricity
    
    #checking if the trajectory is an orbit (ecc<1) or a parabolic/hyperbolic trajectory
    def isOrbit(self):
        if self._eccen < 1.0:
            return True
        else:
            return False
    #semi-major axis (m) if the trajectory is an orbit
    def sem_maj(self):
        if self.isOrbit():
            return -1 * self._grav_pams/(2*self._spec_orb_energ)
    
    #apoapsis radius, in meters
    def apoapsis(self):
        if self.isOrbit():
            return (1 + self._eccen)*self.sem_maj()
    
    #periapsis radius, in meters
    def periapsis(self):
        if self.isOrbit():
            return (1 - self._eccen)*self.sem_maj()
    
    #semi-minor axis (m) if the trajectory is an orbit
    def sem_min(self):
        if self.isOrbit():
            return math.sqrt(self.apoapsis() * self.periapsis())

    #semi-latus rectum (m) if the trajectory is an orbit
    def sem_lat(self):
        if self.isOrbit():
            return self.sem_min()**2 / self.sem_maj()
    
    #orbital period, in seconds
    def period(self):
        if self.isOrbit():
            return 2*math.pi*math.sqrt(self.sem_maj()**3 / self._grav_pams)
    
    #factor by which the simulation is sped up or slowed down compared to real time, all trajectories should be on screen for ten seconds
    def time_factor(self):
        if self.isOrbit():
            return self.period()/10.0

    #mean anomaly, radians
    '''def mean_anom(self, time):
        if self.isOrbit():
            return (2*math.pi/self.period())*time'''

    #eccentric anomaly, radians (uses newtonian approximation based on M = E-esinE
    '''def eccen_anom(self):
        if self.isOrbit():
            M = self.mean_anom()
            ecc = self._eccen
            return (M + ecc*math.sin(M) + (ecc**2)*math.sin(M)*math.cos(M) + 0.5*(ecc**3)*math.sin(M)*(3*(math.cos(M)**2)-1))'''

    #current radius from host, meters
    def radius(self, time):
        if self.isOrbit():
            M = (2*math.pi/self.period())*time
            ecc = self._eccen
            eccen_anom=M + ecc*math.sin(M) + (ecc**2)*math.sin(M)*math.cos(M) + 0.5*(ecc**3)*math.sin(M)*(3*(math.cos(M)**2)-1)
            return(self.sem_maj()*(1-self._eccen*math.cos(eccen_anom)))


#test code, other orbital systems can be entered
earth = Satellite(1.47095e11, 30300, 5.972e24)
sun = Host(1.989e30)
earth_sun=Orbital_System(earth, sun)

#momentary iteration of the system
def sim_iteration(orb_sys, cur_time):
    rad = orb_sys.radius(cur_time*orb_sys.time_factor())
    vel = math.sqrt(orb_sys._grav_pams* ((2/rad) - (1/orb_sys.sem_maj())))

    sys.stdout.write("\rRadius: %f" % (rad/1000) + " kilometers")
    sys.stdout.flush()
    time.sleep(0.000001)   #tells it to print ever 0.1s

#actually running the sim
def run_sim(orb_sys):
    if orb_sys.isOrbit():     #only running if ecc <1, if it's an orbit

        print("Eccentricity: ", orb_sys._eccen)
        print("Semi-major axis: ", orb_sys.sem_maj(), "meters")
        print("Periapsis: ", orb_sys.periapsis(), "meters")
        print("Apoapsis: ", orb_sys.apoapsis(), "meters")
        print("Period: ", orb_sys.period()/86400, "days")
        print("Time acceleration factor: ", orb_sys.time_factor())

        start = time.time()
        while (time.time()-start)>=0.0:
            current_time=time.time()-start
            sim_iteration(orb_sys, current_time)

    #Error message for non-orbits
    else:
        print("Error: Eccentricity >= 1, not an orbit")



run_sim(earth_sun)