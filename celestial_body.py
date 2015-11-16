from math import *
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import *
def _keplerEq(E, eccentricity, mean_anomaly):
    return E + eccentricity*sin(E) - mean_anomaly
def _keplerEqPrime(E, eccentricy, mean_anomaly):
    return 1 + eccentricy*cos(E)
def _keplerEqPrime2(E, eccentricity, mean_anomaly):
    return -eccentricity * sin(E)


class CelestialOrbit:
    def __init__(self, celestial_parent_body, inclination, longitude_ascending_node, argument_periapsis, eccentricity = None, semimajoraxis = None, periapsis = None, apoapsis = None):
        if (eccentricity is not None and semimajoraxis is not None) ^ (periapsis is not None and apoapsis is not None):
            if eccentricity is not None and semimajoraxis is not None:
                self.e = eccentricity
                self.sma = semimajoraxis
            else:
                self.e = (apoapsis - periapsis) / (apoapsis + periapsis)
                self.sma = (apoapsis + periapsis) / 2
            self.parent = celestial_parent_body
            self.i = inclination
            self.long_asc_node = longitude_ascending_node
            self.arg_peri = argument_periapsis
            self.mean_motion = sqrt(self.parent.mu / self.sma**3.)
        else:
            raise TypeError("Incorrect amount of orbital variables")
    def __repr__(self):
        return "CelestialOrbit(celestial_parent_body=%s, inclination=%r, longitude_ascending_node=%r, argument_periapsis = %r, eccentricity=%r, semimajoraxis = %r)" \
            % (self.parent, self.i, self.long_asc_node, self.arg_peri, self.e, self.sma)


    def OrbitalPeriod(self):
        return 2 * pi * sqrt(self.sma ** 3/self.parent.mu)

    def SpecificOrbitalEnergy(self):
        return -self.parent.mu / ( 2 * self.sma)

    def VisViva(self, radius):
        return sqrt(2/radius - 1 / self.sma)

    def RadiusFromTrueAnomaly(self, nu):
        return self.sma * (1. - self.e**2.) / (1 + self.e * cos(nu))
    def RadiusFromEccentricAnomaly(self, E):
        return self.sma * ( 1 - self.e - cos(E))
    def TrueFromEccentricAnomaly(self, E):
        if self.e == 0:
            return E
        return 2 * atan2(sqrt(1+self.e) * sin(E/2), sqrt(1 - self.e) * cos(E/2))
    def EccentricFromTrueAnomaly(self, nu):
        if self.e == 0:
            return nu
        return 2 * atan2(sqrt(1 - self.e) * sin(nu/2), sqrt(1 + self.e) * cos(nu/2))
    def MeanFromEccentricAnomaly(self, E):
        return E - self.e * sin(E)
    def EccentricFromMeanAnomaly(self, M):
        #0 = E + e * sin(E) - M

        guess = M
        is_negative = True
        ecc = self.e
        if ecc == 0:
            return M


        if ecc < .3:
            guess = atan2( sin( M), cos( M) - ecc)
        else:
            if M < 0.:
                M = -M
                is_negative = False

            if ecc > .8 and M < pi / 3. or ecc > 1.:
                trial = M / abs( 1. - ecc)

                if trial * trial > 6. * abs(1. - ecc):
                    if M < pi:
                        trial = ( 6. * M) ** (1./3.)
                    else:
                        trial = asinh( M / ecc)
                guess = trial

        E = newton(_keplerEq, guess, fprime = _keplerEqPrime, fprime2= _keplerEqPrime2 , args=(self.e, M), maxiter = 10)
        return -E if is_negative else E
    def InShadow(self, solar_anomaly):
        avg_theta = (solar_anomaly + pi) %  (2 * pi)
        R_body = self.parent.radius

        def func(theta):
            r = self.Radius(avg_theta + theta)
            return abs(r * sin(theta)) - R_body

        tol = 1.48e-8

        if abs(func(0)) < tol:
            print(0, func(0))
            t2 = avg_theta
            t1 = avg_theta
        elif abs(func(pi/2)) < tol:
            t1 =  avg_theta + pi/2
            t2 = avg_theta + brentq(func, 0, - pi/2, rtol=tol)
        elif abs(func(-pi/2) < tol ):
            t1 = avg_theta + brentq(func, 0, pi/2, rtol=tol)
            t2 = avg_theta -pi/2
        else:
            try:
                t2 = avg_theta + brentq(func, 0, - pi/2, rtol=tol)
                t1 = avg_theta + brentq(func, 0, + pi/2, rtol = tol)
            except ValueError as err:
                print(func(pi/2))
                print("-------------------------------")
                x = np.linspace(-pi, pi, 1E5)
                y=np.array([func(tx) for tx in x])
                plt.plot(x,y)
                plt.show()
                raise

        theta = (t2 , t1 )
        return theta
    def TimeInShadow(self, solar_anomaly):
        theta = self.InShadow(solar_anomaly)

        if theta[1] >= 2*pi:
            theta = (theta[0] - 2*pi, theta[1] - 2*pi)
        E = tuple(self.EccentricFromTrueAnomaly(t) for t in theta)
        M = tuple(self.MeanFromEccentricAnomaly(t) for t in E)

        res = tuple(t/self.mean_motion for t in M)
        #print(theta,E,M, t)
        return res[1]-res[0]

    def MaxTimeInShadow(self):
        if self.e == 0:
            return (0, self.TimeInShadow(0))
        min_obj = minimize_scalar(lambda x: -self.TimeInShadow(x), bounds=(0,pi),method="bounded")
        return min_obj.x, self.TimeInShadow(min_obj.x)

    def TotalMaxTimeInShadow(self):
        orbit = self
        total_time = orbit.MaxTimeInShadow()[1]
        #print("base", total_time,total_time)


        p = orbit.parent
        orbit = orbit.parent.orbit

        while isinstance(orbit.parent, CelestialBody_Orbitting) and orbit.parent.orbit != None:
            t = orbit.MaxTimeInShadow()
            total_time += t[1]
            #print(str(p), t[1],total_time)
            p = orbit.parent
            orbit = orbit.parent.orbit

        return total_time
    def Radius(self, theta):
        return self.sma * (1 - self.e ** 2.) / (1 + self.e * cos(theta))


class CelestialBody:
    def __init__(self, unique_name, mass, mu, radius, sidereal_period, atmosphere_height = 0, low_orbit_height = None):
        self.mass = mass
        self.mu = mu
        self.radius = radius
        self.sidereal_period = sidereal_period
        self.name = unique_name
        self.childs = set()
        self.atmosphere_height = atmosphere_height
        if low_orbit_height is None:
            if self.atmosphere_height == 0:
                self.low_orbit_height = self.radius + 20000
            else:
                self.low_orbit_height = self.radius + self.atmosphere_height + 10000
        else:
            self.low_orbit_height = low_orbit_height
    def __repr__(self):
        return "CelestialBody(unique_name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r)" \
            % (self.name,self.mass, self.mu, self.radius, self.sidereal_period, self.atmosphere_height, self.low_orbit_height)
    def __str__(self):
        return self.name


    def AppendChild(self, celestial_body):
        self.childs.add(celestial_body)

    def DiscardChild(self, celestial_body):
        self.childs(celestial_body)

    def SurfaceSpeed(self, longitude):
        r = self.radius * cos(radians(longitude))
        l = 2 * pi * r

        v = l / self.radius
        return v

    def CreateGeoStatOrbit(self):
        i = 0
        long_asc = 0
        arg_peri = 0
        #T = = 2 * pi * sqrt(a ^ 3 / mu)
        #a = ((T/(2*pi))^2 * mu) ^ (1/3)
        sma = ((self.sidereal_period / (2 * pi)) ** 2. * self.mu)
        return CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)

    def CreateLowOrbit(self):
        sma = self.radius + self.low_orbit_height
        i = 0
        long_asc = 0
        arg_peri = 0
        return CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)








class CelestialBody_Orbitting(CelestialBody):
    def __init__(self, unique_name, mass, mu, radius, sidereal_period, SOI_radius, atmosphere_height = 0, low_orbit_height = None, mean_anomaly_at_epoch = 0, orbit = None):
        super().__init__(unique_name, mass, mu, radius, sidereal_period, atmosphere_height, low_orbit_height)
        if orbit is None:
            self.orbit = None
        else:
            self.orbit = orbit
            p = self.orbit.parent
            p.AppendChild(self)
        self.SOI = SOI_radius
        self.M0 = mean_anomaly_at_epoch
    def __repr__(self):
        return "CelestialBody_Orbitting(unique_name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r, mean_anomaly_at_epoch = %r, \n\t orbit = %r)" \
            % (self.name,self.mass, self.mu, self.radius, self.sidereal_period, self.atmosphere_height, self.low_orbit_height, self.M0, self.orbit)


    def getMeanAnomaly(self, time):
        return self.M0 + self.orbit.mean_motion * time

    def getEccentricAnomaly(self, time):
        M = self.getMeanAnomaly(time)
        return self.orbit.EccentricFromMeanAnomaly(M)

    def getTrueAnomaly(self, time):

        M = self.getMeanAnomaly(time)
        E = self.orbit.EccentricFromMeanAnomaly(M)
        return self.orbit.TrueFromEccentricAnomaly(E)

    def RegisterOrbit(self, orbit):
        if self.orbit is not None:
            p = self.orbit.parent
            p.DiscardChild(self)
        self.orbit = orbit
        p = self.orbit.parent
        p.AppendChild(self)






