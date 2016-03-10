from math import pi
import numpy as np
import math as m
from numpy import cos, sin, sqrt, power, square, arctan2, arccos, arcsin, arcsinh, radians, degrees
import celestial_orbit as CO
import typing

atan2 = arctan2
acos = arccos
asin = arcsin
asinh = arcsinh


class CelestialObjectAbstract:
    def __init__(self, name:str, current_time:float=0.0, pos:np.ndarray=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name  # type:str
        self.current_time = current_time  # type:foat
        if pos is None:
            self.__local_pos = np.array([0, 0, 0])  # type:np.ndarray
        else:
            self.__local_pos = np.array(pos)  # type:np.ndarray
        self.__last_calculated_time = None  # type:float
        self.__global_pos = np.array([0, 0, 0])  # type:np.ndarray

    def __repr__(self):
        return "CelestialBody_NoProperty(name=%r)" \
               % (self.name)

    def __str__(self):
        return self.name

    def getRoot(self) -> "CelestialObjectAbstract":
        o = self
        try:
            while True:
                o = self.parent
        except AttributeError:
            return o

    def _updateTimeWholeSystem(self, new_time:float, calculate_now:bool=False):
        self.current_time = new_time
        if calculate_now:
            self.get_global_position_at_time(new_time, True)
        for body in self.childs:
            body._updateTimeWholeSystem(new_time, calculate_now)

    def updateTimeWholeSystem(self, new_time:float, calculate_now:bool=False):
        # get root
        root = self.getRoot()
        root._updateTimeWholeSystem(new_time, calculate_now)

    def __MakePositionDirty(self):
        self.__last_calculated_time = None
        for c in self.childs:
            c.__MakePositionDirty()

    def set_local_position(self, newpos:np.ndarray):
        self.__local_pos = newpos

    def get_local_position_at_time(self, time:float=None):
        return self.__local_pos

    def get_global_position_at_time(self, time:float=None, do_update:bool=False) -> np.ndarray:
        if time is not None and self.__last_calculated_time == time:
            return self.__global_pos
        if time is None:
            time = self.current_time
        pos = self.localToGlobalPosition(self.get_local_position_at_time(time), time, do_update)
        if do_update:
            self.current_time = time
            self.__global_pos = pos
        return pos

    def localToGlobalPosition(self, position, time=None, do_update=False):
        try:
            obj = self.orbit.parent
            par_pos = obj.get_global_position_at_time(time, do_update)
            return position + par_pos
        except AttributeError:
            return position

    def clear_cached_location(self):
        self.__last_calculated_time = None


class CelestialObjectAbstract_Orbitting(CelestialObjectAbstract):
    def __init__(self, mean_anomaly_at_epoch: CO.MeanAnomaly = 0., orbit: CO.CelestialOrbit = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._orbit = orbit  # type:CO.CelestialOrbit
        if orbit is not None:
            p = self.parent
            p.AppendChild(self)
        self.M0 = CO.MeanAnomaly(mean_anomaly_at_epoch)  # type:CO.MeanAnomaly
        self.__last_calculated_mean_anomaly = None  # type:CO.MeanAnomaly
        self.__stored_eccentric_anomaly = None  # type:CO.EccentricAnomaly

    def get_meananomaly(self, time: float = None) -> CO.MeanAnomaly:
        if time is None:
            time = self.current_time
        return self.M0 + CO.MeanAnomaly(self.orbit.mean_motion * time)

    def get_eccentricanomaly(self, time:float=None) -> CO.EccentricAnomaly:
        if time is None:
            time = self.current_time
        M = self.get_meananomaly(time)
        if M != self.__last_calculated_mean_anomaly:
            E = self.orbit.eccentric_from_mean_anomaly(M)
            self.__last_calculated_mean_anomaly = M
            self.__stored_eccentric_anomaly = E
        else:
            E = self.__stored_eccentric_anomaly
        return E

    def get_trueanomaly(self, time:float=None) -> CO.TrueAnomaly:
        E = self.get_eccentricanomaly(time)
        return self.orbit.true_from_eccentric_anomaly(E)

    def register_orbit(self, orbit:CO.CelestialOrbit):
        if self.orbit is not None:
            p = self.parent
            p.DiscardChild(self)
        self._orbit = orbit
        self.__last_calculated_mean_anomaly = None
        if orbit is not None:
            p = self.parent
            p.AppendChild(self)

    @property
    def parent(self) -> "CelestialBody":
        return self.orbit.parent

    def get_distance_central_body(self, time:float=None) -> float:
        E = self.get_eccentricanomaly(time)
        d = self.orbit.radius_from_eccentric_anomaly(E)
        return d

    def get_local_postion_at_trueanomaly(self, true_anomaly: CO.TrueAnomaly) -> np.ndarray:
        return self.orbit.true_anomaly_to_position(true_anomaly)

    def get_local_position_at_time(self, time: float = None) -> np.ndarray:
        if time is None:
            time = self.__current_time
        return self.get_local_postion_at_trueanomaly(self.get_trueanomaly(time))

    @property
    def orbit_period(self) -> float:
        return self.orbit.orbital_period()

    @property
    def orbit(self) -> CO.CelestialOrbit:
        return self._orbit


class CelestialBody(CelestialObjectAbstract):
    def __init__(self, mass: float, mu: float, radius: float, sidereal_rotation_period: float,
                 SOI_radius: float = m.inf,
                 atmosphere_height: float = 0., low_orbit_height: float = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mass = mass  # type: float
        self.mu = mu  # type: float
        self.radius = radius  # type:float
        self.sidereal_rotation_period = sidereal_rotation_period  # type:float
        self.childs = set()  # type:typing.Set
        self.atmosphere_height = atmosphere_height  # type:float
        self.SOI = SOI_radius  # type:float

        if low_orbit_height is None:
            if self.atmosphere_height == 0:
                self.low_orbit_height = self.radius + 25000  # type:float
            else:
                self.low_orbit_height = self.radius + self.atmosphere_height + 10000  # type:float
        else:
            self.low_orbit_height = low_orbit_height  # type:float

    def __repr__(self):
        return "CelestialBody(name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r)" \
               % (self.name, self.mass, self.mu, self.radius, self.sidereal_rotation_period, self.atmosphere_height,
                  self.low_orbit_height)

    def __str__(self):
        return self.name

    def AppendChild(self, celestial_body:CelestialObjectAbstract):
        """
        Appends a satellite to the parent
        :param celestial_body: child
        :rtype: None
        """
        self.childs.add(celestial_body)

    def DiscardChild(self, celestial_body:CelestialObjectAbstract):
        self.childs.remove(celestial_body)

    def SurfaceSpeed(self, longitude:float) -> float:
        r = self.radius * cos(radians(longitude))
        l = 2 * pi * r
        v = l / self.sidereal_rotation_period
        return v

    def CreateGeoStatOrbit(self) -> CO.CelestialOrbit:
        i = 0
        long_asc = 0
        arg_peri = 0
        # T = = 2 * pi * sqrt(a ^ 3 / mu)
        # a = ((T/(2*pi))^2 * mu) ^ (1/3)
        sma = ((self.sidereal_rotation_period / (2 * pi)) ** 2. * self.mu) ** (1. / 3.)
        return CO.CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)

    def CreateLowOrbit(self) -> CO.CelestialOrbit:
        sma = self.radius + self.low_orbit_height
        i = 0
        long_asc = 0
        arg_peri = 0
        return CO.CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)


class CelestialStar(CelestialBody):
    def __init__(self, brightness:float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__brightness = brightness #type:float

    @property
    def brightness(self) -> float:
        return self.__brightness

    @brightness.setter
    def brightness(self, value:float):
        self.__brightness = value

    def getPower(self, distance:float) -> float:
        return self.__brightness / distance ** 2


class CelestialBody_Orbitting(CelestialBody, CelestialObjectAbstract_Orbitting):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return "CelestialBody_Orbitting(unique_name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r, mean_anomaly_at_epoch = %r, \n\t orbit = %r)" \
               % (self.name, self.mass, self.mu, self.radius, self.sidereal_rotation_period, self.atmosphere_height,
                  self.low_orbit_height, self.M0, self.orbit)


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
