from math import pi
import numpy as np
from numpy import cos, sin, sqrt, power, square, arctan2, arccos, arcsin, arcsinh, radians, degrees
import matplotlib.pyplot as plt
from scipy.optimize import *
import scipy as sp
import celestial_orbit as CO

atan2 = arctan2
acos = arccos
asin = arcsin
asinh = arcsinh




class CelestialBody(object):
    def __init__(self, name, mass, mu, radius, sidereal_rotation_period, atmosphere_height=0., low_orbit_height=None, pos = None, current_time = 0.0, *args, **kwargs):
        """
        :type name: str
        :type mass: float
        :type mu: float
        :type radius: float
        :type sidereal_rotation_period: float
        :type atmosphere_height: float
        :type low_orbit_height: float
        :type pos: np.ndarray
        :type current_time: float
        :param name: reference
        :param mass: mass
        :param mu: gravitational parameter
        :param radius: radius
        :param sidereal_rotation_period: rotation
        :param atmosphere_height: height atmosphere
        :param low_orbit_height: height low orbit
        :param pos: position in space (default = 0,0,0)
        :param current_time: time for calculations
        :return:
        """
        super().__init__(*args, **kwargs)
        self.mass = mass
        self.mu = mu
        self.radius = radius
        self.sidereal_rotation_period = sidereal_rotation_period
        self.name = name
        self.childs = set()
        self.atmosphere_height = atmosphere_height
        self.current_time = current_time

        if low_orbit_height is None:
            if self.atmosphere_height == 0:
                self.low_orbit_height = self.radius + 20000
            else:
                self.low_orbit_height = self.radius + self.atmosphere_height + 10000
        else:
            self.low_orbit_height = low_orbit_height
        if pos is None:
            self.__local_pos = np.array([0,0,0])
        else:
            self.__local_pos = np.array(pos)

        self.__last_calculated_time = None
        self.__global_pos = np.array([0,0,0])

    def __repr__(self):
        return "CelestialBody(unique_name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r)" \
               % (self.name, self.mass, self.mu, self.radius, self.sidereal_rotation_period, self.atmosphere_height,
                  self.low_orbit_height)

    def __str__(self):
        return self.name

    def AppendChild(self, celestial_body):
        """
        Appends a satellite to the parent
        :param celestial_body: child
        :type celestial_body: CelestialBody
        :rtype: None
        """
        self.childs.add(celestial_body)

    def DiscardChild(self, celestial_body):
        self.childs.remove(celestial_body)

    def SurfaceSpeed(self, longitude):
        r = self.radius * cos(radians(longitude))
        l = 2 * pi * r

        v = l / self.radius
        return v

    def CreateGeoStatOrbit(self):
        i = 0
        long_asc = 0
        arg_peri = 0
        # T = = 2 * pi * sqrt(a ^ 3 / mu)
        # a = ((T/(2*pi))^2 * mu) ^ (1/3)
        sma = ((self.sidereal_rotation_period / (2 * pi)) ** 2. * self.mu) ** (1. / 3.)
        return CO.CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)

    def CreateLowOrbit(self):
        sma = self.radius + self.low_orbit_height
        i = 0
        long_asc = 0
        arg_peri = 0
        return CO.CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)

    def getRoot(self):
        o = self
        try:
            while True:
                o = self.orbit.parent
        except AttributeError:
            return o

    def _updateTimeWholeSystem(self, new_time, calculate_now = False):
        self.current_time = new_time
        if calculate_now:
            self.getGlobalPositionAtTime(new_time, True)
        for body in self.childs:
            body._updateTimeWholeSystem(new_time, calculate_now)


    def updateTimeWholeSystem(self, new_time, calculate_now = False):
        #get root
        root = self.getRoot()
        root._updateTimeWholeSystem(new_time, calculate_now)



    def __MakePositionDirty(self):
        self.__last_calculated_time = None
        for c in self.childs:
            c.__MakePositionDirty()

    def setLocalPosition(self, newpos):
        self.__local_pos = newpos

    def getLocalPositionAtTime(self, time=None):
        return self.__local_pos

    def getGlobalPositionAtTime(self, time=None, do_update = False):
        if time is not None and self.__last_calculated_time == time:
            return self.__global_pos
        if time is None:
            time = self.current_time
        pos = self.localToGlobalPosition(self.getLocalPositionAtTime(time), time, do_update)
        if do_update:
            self.current_time = time
            self.__global_pos = pos
        return pos

    def localToGlobalPosition(self, position, time = None, do_update = False):
        try:
            obj = self.orbit.parent
            par_pos = obj.getGlobalPositionAtTime(time, do_update)
            return position + par_pos
        except AttributeError:
            return position

    def clearCachedLocation(self):
        self.__last_calculated_time = None

class CelestialStar(CelestialBody):
    def __init__(self, brightness, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__brightness = brightness
    @property
    def brightness(self):
        return self.__brightness

    @brightness.setter
    def brightness(self, value):
        self.__brightness = value

    def getPower(self, distance):
        self.__brightness/distance**2


class CelestialBody_Orbitting(CelestialBody):
    def __init__(self, SOI_radius, mean_anomaly_at_epoch=0, orbit=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if orbit is None:
            self._orbit = None
        else:
            self._orbit = orbit
            p = self.orbit.parent
            p.AppendChild(self)
        self.SOI = SOI_radius
        self.M0 = CO.MeanAnomaly(mean_anomaly_at_epoch)
        self.__last_calculated_mean_anomaly = None
        self.__stored_eccentric_anomaly = None
    def __repr__(self):
        return "CelestialBody_Orbitting(unique_name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r, mean_anomaly_at_epoch = %r, \n\t orbit = %r)" \
               % (self.name, self.mass, self.mu, self.radius, self.sidereal_rotation_period, self.atmosphere_height,
                  self.low_orbit_height, self.M0, self.orbit)

    def getMeanAnomaly(self, time = None):
        if time is None:
            time = self.current_time
        return self.M0 + CO.MeanAnomaly(self.orbit.mean_motion * time)

    def getEccentricAnomaly(self, time = None):
        if time is None:
            time = self.current_time
        M = self.getMeanAnomaly(time)
        if M != self.__last_calculated_mean_anomaly:
            E = self.orbit.eccentric_from_mean_anomaly(M)
            self.__last_calculated_mean_anomaly = M
            self.__stored_eccentric_anomaly = E
        else:
            E = self.__stored_eccentric_anomaly
        return E

    def getTrueAnomaly(self, time = None):
        E = self.getEccentricAnomaly(time)
        return self.orbit.true_from_eccentric_anomaly(E)

    def RegisterOrbit(self, orbit):
        if self.orbit is not None:
            p = self.orbit.parent
            p.DiscardChild(self)
        self.orbit = orbit
        p = self.orbit.parent
        p.AppendChild(self)

    def get_distance_central_body(self, time = None):
        E = self.getEccentricAnomaly(time)
        d = self.orbit.radius_from_eccentric_anomaly(E)
        return d


    def getLocalPositionAtTrueAnomaly(self, true_anomaly):
        return self.orbit.true_anomaly_to_position(true_anomaly)
    def getLocalPositionAtTime(self, time=None):
        if time is None:
            time = self.__current_time
        return self.getLocalPositionAtTrueAnomaly(self.getTrueAnomaly(time))

    @property
    def orbit_period(self):
        return self.orbit.orbital_period()
    @property
    def orbit(self):
        return self._orbit

    @orbit.setter
    def orbit(self, value):
        self.orbit = value
        self.__last_calculated_mean_anomaly = None

class CelestialBody_OrbittingStar(CelestialBody_Orbitting, CelestialStar):
    pass



def plotSphere(x, y, z, r, ax, steps=20, **kwargs):
    # u, v = np.mgrid[0:2*pi:20j, 0:np.pi:10j]
    # x+=r*cos(u)*sin(v)
    # y+=r*sin(u)*sin(v)
    # z+=r*cos(v)    #ax.plot_wireframe(x, y, z, **kwargs)

    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(0, np.pi, steps)

    x += r * np.outer(np.cos(u), np.sin(v))
    y += r * np.outer(np.sin(u), np.sin(v))
    z += r * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, rstride=1, cstride=1, **kwargs)
