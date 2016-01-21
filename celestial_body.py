from math import pi
import numpy as np
from numpy import cos, sin, sqrt, power, square, arctan2, arccos, arcsin, arcsinh, radians, degrees
import matplotlib.pyplot as plt
import sys
from scipy.optimize import *

atan2 = arctan2
acos = arccos
asin = arcsin
asinh = arcsinh


def _keplerEq(E, eccentricity, mean_anomaly):
    return E - eccentricity * sin(E) - mean_anomaly


def _keplerEqPrime(E, eccentricy, mean_anomaly):
    return 1 - eccentricy * cos(E)


def _keplerEqPrime2(E, eccentricity, mean_anomaly):
    return eccentricity * sin(E)


def _TransformCoordinates(mat, vec):
    return mat.dot(vec)


def _VectorRotateOriginAxis(vec, unit_dir_vector, theta):
    u, v, w = unit_dir_vector
    x, y, z = vec
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    P = (u * x + v * y + w * z) * (1 - cos_theta)
    res = np.array([
        u * P + x * cos_theta + (-w * y + v * z) * sin_theta,
        v * P + y * cos_theta + (w * x - u * z) * sin_theta,
        w * P + z * cos_theta + (-v * x + u * y) * sin_theta
    ])
    return res


def _VectorRotateX(vec, theta):
    mat = np.array([[1, 0, 0],
                    [0, cos(theta), -sin(theta)],
                    [0, sin(theta), cos(theta)]])
    return mat.dot(vec)


def _VectorRotateY(vec, theta):
    mat = np.array([[cos(theta), 0, -sin(theta)],
                    [0, 1, 0],
                    [sin(theta), 0, cos(theta)]])
    return mat.dot(vec)


def _VectorRotateZ(vec, theta):
    mat = np.array([[cos(theta), -sin(theta), 0],
                    [sin(theta), cos(theta), 0],
                    [0, 0, 1]])
    print(mat)
    return mat.dot(vec)


class Anomaly(float):
    """
    Anomaly class
    """

    @staticmethod
    def FromDegrees(deg):
        """
        :param deg: degrees
        :return: anomaly in radians
        :rtype: Anomaly
        """
        return Anomaly(radians(deg))

    def ToDegrees(self):
        return degrees(self)

    def ModulateFullCircle(self):
        """
        Modulates value to [0, 2*pi) interval
        :return: Modulated anomaly
        """
        return type(self)(self % (2 * pi))

    def ModulateBidirectional(self):
        """
        Modulates value to (-pi, pi] interval
        :return: Modulated anomaly
        """
        p = self % (2 * pi)
        if p > pi:
            p -= 2 * pi
        type(self)(p)


class TrueAnomaly(Anomaly):
    """
    True Anomaly class
    """


class EccentricAnomaly(Anomaly):
    """
    True Anomaly class
    """


class MeanAnomaly(Anomaly):
    """
    Mean anomaly class
    """


class CelestialOrbit:
    """
    An orbit around a given celestial body. Multiple bodies can use the same orbit.
    It contains all elements apart from the anomaly

    Attributes
    ----------
    e: eccentricity
    sma: semi major axis
    parent: central body
    i: inclination
    long_asc_node: longitude of ascending node
    arg_per: argument of periapsis
    mean_motion: mean motion
    """

    def __init__(self, celestial_parent_body, inclination, longitude_ascending_node, argument_periapsis,
                 eccentricity=None, semimajoraxis=None, semilatusrectum=None, periapsis=None, apoapsis=None):
        """
        Creates an orbit around a given celestial body. Multiple bodies can use the same orbit.
        It contains all elements apart from the anomaly
        Either provide both eccentricity & semimajoraxis or periapsis & apoapsis
        The two that are not given are calculated from the others
        :type celestial_parent_body: CelestialBody
        :type inclination: float
        :type longitude_ascending_node: float
        :type argument_periapsis: float
        :type eccentricity: float
        :type semimajoraxis: float
        :type semilatusrectum: float
        :type periapsis: float
        :type apoapsis: float

        :param celestial_parent_body: Parent body
        :param inclination: Inclination
        :param longitude_ascending_node: Longitude of ascending node
        :param argument_periapsis: Argument of periapsis
        :param eccentricity: Eccentricity
        :param semimajoraxis: Semimajor axis
        :param semilatusrectum: semilatus rectum
        :param periapsis: Periapsis
        :param apoapsis: Apoapsis
        :return: Celestial orbit
        """
        self.setOrbitSize(eccentricity, semimajoraxis, semilatusrectum, periapsis, apoapsis)
        self.parent = celestial_parent_body
        self.__i = inclination
        self.__long_asc_node = longitude_ascending_node
        self.__arg_peri = argument_periapsis
        self.__UpdateConvMats()

    def __repr__(self):
        return "CelestialOrbit(celestial_parent_body=%s, inclination=%r, longitude_ascending_node=%r, argument_periapsis = %r, eccentricity=%r, semimajoraxis = %r)" \
               % (self.parent, self.inclination, self.long_asc_node, self.arg_peri, self.e, self.sma)

    def OrbitalPeriod(self):
        """
        T = 2*pi * sqrt(a^3 / mu)
        :rtype: float
        :return: Period
        :rtype: float
        """
        return 2 * pi * sqrt(self.sma ** 3 / self.parent.mu)

    def SpecificOrbitalEnergy(self):
        """
        E = -mu/(2*a)
        :return: Specific Orbital energy [Joule]
        :rtype: float
        """
        return -self.parent.mu / (2 * self.sma)

    def VisViva(self, radius):
        """
        Returns speed given a distance

        :param radius: Radius
        :type radius: float
        :rtype: float
        :return: speed
        """
        return sqrt(2 / radius - 1 / self.sma)

    def RadiusFromTrueAnomaly(self, nu):
        """
        :param nu: True anomaly
        :type nu: TrueAnomaly
        :return: Radius
        :rtype: float
        """
        return self.__semi_latus_rectum / (1 + self.__e * cos(nu))

    def RadiusFromEccentricAnomaly(self, E):
        """
        :param E: Eccentric anomaly
        :type E: EccentricAnomaly
        :return: Radius
        :rtype: float
        """
        return self.sma * (1 - self.e - cos(E))

    def ToTrueAnomaly(self, theta):
        """
        :param theta: Anomaly
        :type theta: Anomaly
        :return: True Anomaly
        :rtype: TrueAnomaly
        """
        if isinstance(theta, EccentricAnomaly):
            return self.TrueFromEccentricAnomaly(theta)
        elif isinstance(theta, MeanAnomaly):
            return self.TrueFromEccentricAnomaly(self.EccentricFromMeanAnomaly(theta))
        else:
            return TrueAnomaly(theta)

    def TrueFromEccentricAnomaly(self, E):
        """
        :param E: Eccentric Anomaly
        :type E: EccentricAnomaly
        :return: True Anomaly
        :rtype: TrueAnomaly
        """
        if self.e == 0:
            return E
        try:
            return TrueAnomaly(2 * atan2(sqrt(1 + self.e) * sin(E / 2), sqrt(1 - self.e) * cos(E / 2)))
        except TypeError:
            return np.array(2 * atan2(sqrt(1 + self.e) * sin(E / 2), sqrt(1 - self.e) * cos(E / 2)))

    def ToEccentricAnomaly(self, theta):
        """
        :param theta: Anomaly
        :type theta: Anomaly
        :return: Eccentric Anomaly
        :rtype: EccentricAnomaly
        """
        if isinstance(theta, TrueAnomaly):
            return self.EccentricFromTrueAnomaly(theta)
        elif isinstance(theta, MeanAnomaly):
            return self.EccentricFromMeanAnomaly(theta)
        else:
            return EccentricAnomaly(theta)

    def EccentricFromTrueAnomaly(self, nu):
        """
        :type nu: TrueAnomaly
        :param nu: True anomaly
        :return: Eccentric Anomaly
        :rtype: EccentricAnomaly
        """
        if self.e == 0:
            return nu
        return 2 * atan2(sqrt(1 - self.e) * sin(nu / 2), sqrt(1 + self.e) * cos(nu / 2))

    def _EccentricFromMeanAnomalyNewtonMethod(self, M):
        """
        Gets the Eccentric anomaly from mean anomaly using Newton's numerical root finding
        :type M: MeanAnomaly
        :param M: Mean anomaly
        :return: Eccentric Anomaly
        :rtype: EccentricAnomaly
        """

        # 0 = E + e * sin(E) - M
        is_negative = False
        ecc = self.e
        if ecc == 0:
            return M

        M %= 2 * pi
        if M >= pi:
            M = M - 2 * pi
        elif M < -pi:
            M = 2 * pi + M
        if ecc < .3:
            guess = atan2(sin(M), cos(M) - ecc)
        else:
            if M < 0.:
                M = -M
                is_negative = True
            guess = M
            if ecc > .8 and M < pi / 3. or ecc > 1.:
                trial = M / abs(1. - ecc)

                if trial * trial > 6. * abs(1. - ecc):
                    if M < pi:
                        trial = (6. * M) ** (1. / 3.)
                    else:
                        trial = asinh(M / ecc)
                guess = trial

        try:
            # E = newton(_keplerEq, guess, fprime = _keplerEqPrime, fprime2= _keplerEqPrime2 , args=(self.e, M), maxiter = 100)
            E = newton(_keplerEq, guess, fprime=_keplerEqPrime, args=(self.e, M), maxiter=1000000)
        except RuntimeError as e:
            print("Failed")
            print(repr(self))
            print(-M if is_negative else M, e)
            raise

        return -E if is_negative else E

    def EccentricFromMeanAnomaly(self, M):
        """
        Gets the Eccentric anomaly from mean anomaly using Brent's numerical root finding
        :type M: MeanAnomaly
        :param M: Mean anomaly
        :return: Eccentric Anomaly
        :rtype: EccentricAnomaly
        """
        try:
            n = M // pi

            if M == 0: return 0
            if M == pi: return pi
            E = brentq(_keplerEq, n * pi, (n + 1) * pi, args=(self.e, M), maxiter=1000, disp=True,
                       xtol=np.finfo(float).eps)
            # E = newton(_keplerEq, guess, fprime = _keplerEqPrime, args=(self.e, M), maxiter = 100)
        except RuntimeError as e:
            print("Failed")
            print(repr(self))
            print(M, e)
            raise
        except ValueError as e:
            print(M, n * pi, (n + 1) * pi)
            print(_keplerEq(0, self.e, M), _keplerEq(-pi, self.e, M), _keplerEq(pi, self.e, M))
            print(repr(self))
            raise
        return E

    def ToMeanAnomaly(self, theta):
        """
        :type theta: Anomaly
        :param theta: Anomaly
        :return: Mean Anomaly
        :rtype: MeanAnomaly
        """
        if isinstance(theta, EccentricAnomaly):
            return self.MeanFromEccentricAnomaly(theta)
        elif isinstance(theta, TrueAnomaly):
            return self.MeanFromEccentricAnomaly(self.EccentricFromTrueAnomaly(theta))
        else:
            return MeanAnomaly(theta)

    def MeanFromEccentricAnomaly(self, E):
        """
        :type E: EccentricAnomaly
        :param E: Eccentric anomaly
        :return: Mean Anomaly
        :rtype: MeanAnomaly
        """
        return E - self.e * sin(E)

    def InShadow(self, solar_anomaly):
        """
        :param solar_anomaly: Solar direction
        :type solar_anomaly: Anomaly
        :rtype: (TrueAnomaly, TrueAnomaly)
        :return: Start and ending true anomly of the shadow
        """
        avg_theta = (solar_anomaly + pi) % (2 * pi)
        R_body = self.parent.radius

        def shadow_func(theta):
            r = self.RadiusFromTrueAnomaly(avg_theta + theta)
            return abs(r * sin(theta)) - R_body

        tol = 1.48e-9
        try:
            if abs(shadow_func(0)) < tol:
                t2 = avg_theta
                t1 = avg_theta
            elif abs(shadow_func(pi / 2)) < tol:
                t1 = avg_theta + pi / 2
                t2 = avg_theta + brentq(shadow_func, 0, - pi / 2)
            elif abs(shadow_func(-pi / 2) < tol):
                t1 = avg_theta + brentq(shadow_func, 0, pi / 2)
                t2 = avg_theta - pi / 2
            else:
                t2 = avg_theta + brentq(shadow_func, 0, - pi / 2)
                t1 = avg_theta + brentq(shadow_func, 0, + pi / 2)
        except ValueError as err:
            print("input: ", avg_theta)
            print("bounds: ", shadow_func(-pi / 2), shadow_func(0), shadow_func(pi / 2))
            raise

        theta = (t2, t1)
        return theta

    def TimeInShadow(self, solar_anomaly):
        """
        :param solar_anomaly: Solar direction
        :type solar_anomaly: Anomaly
        :rtype: float
        :return: Start and ending true anomly of the shadow
        """
        theta = self.InShadow(solar_anomaly)

        if theta[1] >= 2 * pi:
            theta = (theta[0] - 2 * pi, theta[1] - 2 * pi)
        E = tuple(self.EccentricFromTrueAnomaly(t) for t in theta)
        M = tuple(self.MeanFromEccentricAnomaly(t) for t in E)

        res = tuple(t / self.mean_motion for t in M)
        return res[1] - res[0]

    def MaxTimeInShadow(self):
        """
        Returns maximum time in shadow.
        Numerical solution, for standard orbits TimeInShadow(0) might be faster
        :return: time
        :rtype: float
        """
        #if self.e == 0:
        return self.TimeInShadow(0)
        #min_obj = minimize_scalar(lambda x: -self.TimeInShadow(x), bounds=(0, pi), method="bounded")
        #return min_obj.x, self.TimeInShadow(min_obj.x)

    def TotalMaxTimeInShadow(self):
        """
        Returns total max time in shadow
        Iterativelly calls MaxTimeInShadow for each orbit until star is found
        :return: time
        :rtype: float
        """
        return sum( i[1] for i in self.GenerateShadowLightTimeList())

    def GenerateShadowLightTimeList(self):
        """
        Returns total max time in shadow
        Iterativelly calls MaxTimeInShadow for each orbit until star is found
        :return: list of (light, shadow) times
        :rtype: [tuple(float, float)]
        """
        all_times=[]
        o = self
        try:
            while not isinstance(o.parent, CelestialStar):
                period = o.period
                t = o.MaxTimeInShadow()
                yield (period-t, t)
                o = o.parent.orbit
        except AttributeError:
            pass
        return all_times

    def PythagoralDistance(self, nu1, nu2):
        """
        :type nu1: TrueAnomaly
        :param nu1: Anomaly first object
        :param nu2: Anomaly second object
        :type nu2: TrueAnomaly
        :return: Distance between objects
        :rtype: float
        """
        p1 = self.TrueAnomalyToPositionVector(nu1)
        p2 = self.TrueAnomalyToPositionVector(nu2)
        return np.linalg.norm(p1 - p2)

    def TrueAnomalyToPositionVector(self, nu):
        """
        :param nu: TrueAnomaly
        :type nu: TrueAnomaly
        :return: Vector
        :rtype: np.multiarray.ndarray
        """
        r_abs = self.RadiusFromTrueAnomaly(nu)
        r_rotated = np.array([r_abs * cos(nu), r_abs * sin(nu), 0, 1])
        r = self.__planar_to_local_mat.dot(r_rotated)
        r = np.squeeze(np.vstack(r))
        return r[:3]

    def TrueAnomalyToVelocityVector(self, nu):
        v_inplane = sqrt(self.parent.mu/self.__semi_latus_rectum) * np.array([-sin(nu), self.e + cos(nu), 0, 0])
        v = self.__planar_to_local_mat.dot(v_inplane)
        return np.squeeze(np.vstack(v))

    def setOrbitSize(self, eccentricity=None, semimajoraxis=None, semilatusrectum=None, periapsis=None, apoapsis=None):
        t = [eccentricity is not None, semimajoraxis is not None, semilatusrectum is not None, periapsis is not None,
             apoapsis is not None]
        str = "Incorrect amount of orbital variables"
        if sum(t) != 2:
            raise TypeError(str)
        if eccentricity is not None:
            self.__e = eccentricity
            if semilatusrectum is not None:
                self.__semi_latus_rectum = semilatusrectum
            elif semimajoraxis is not None:
                self.sma = semimajoraxis
            elif apoapsis is not None:
                self.sma = semimajoraxis
            elif periapsis is not None:
                self.__semi_latus_rectum = periapsis * (eccentricity + 1)
            else:
                raise TypeError(str)
        elif semilatusrectum is not None:
            if periapsis is not None:
                self.__e = semilatusrectum / periapsis - 1
                raise TypeError(str)
        elif semimajoraxis is not None:
            self.sma = semimajoraxis
            if apoapsis is not None:
                periapsis = 2 * semimajoraxis - apoapsis
            if periapsis is not None:
                self.__e = 1 - periapsis / semimajoraxis
            else:
                raise TypeError(str)

    @property
    def inclination(self):
        return self.__i

    @property
    def long_asc_node(self):
        return self.__long_asc_node

    @property
    def arg_peri(self):
        return self.__arg_peri

    @property
    def sma(self):
        return self.__semi_latus_rectum / (1 - self.__e ** 2)

    @sma.setter
    def sma(self, value):
        self.semi_latus_rectum = value * (1 - self.e ** 2)

    @property
    def semi_latus_rectum(self):
        return self.__semi_latus_rectum

    @semi_latus_rectum.setter
    def semi_latus_rectum(self, val):
        self.__semi_latus_rectum = val

    @property
    def e(self):
        return self.__e

    @e.setter
    def e(self, val):
        self.__e = val

    @property
    def mean_motion(self):
        return sqrt(self.parent.mu / self.sma ** 3.)

    @property
    def period(self):
        return 2 * pi * sqrt(self.sma ** 3 / self.parent.mu)

    @property
    def periapsis_distance(self):
        return self.__semi_latus_rectum / (1+self.e)

    @property
    def apoapsis_distance(self):
        return self.__semi_latus_rectum / (1-self.__e)

    def __UpdateConvMats(self):
        t = np.array([[cos(self.long_asc_node), -sin(self.long_asc_node), 0, 0],
                      [sin(self.long_asc_node), cos(self.long_asc_node), 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
        t = t.dot([[1, 0, 0, 0],
                   [0, cos(self.inclination), -sin(self.inclination),0],
                   [0, sin(self.inclination), cos(self.inclination),0],
                   [0, 0, 0, 1]])
        t = t.dot([[cos(self.long_asc_node), -sin(self.long_asc_node), 0, 0],
                      [sin(self.long_asc_node), cos(self.long_asc_node), 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
        self.__planar_to_local_mat = t
        self.__local_to_planar_mat = t.T

class CelestialBody(object):
    def __init__(self, name, mass, mu, radius, sidereal_rotation_period, atmosphere_height=0., low_orbit_height=None, pos = None, current_time = 0.0, **kwargs):
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
        super().__init__(**kwargs)
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
        return CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)

    def CreateLowOrbit(self):
        sma = self.radius + self.low_orbit_height
        i = 0
        long_asc = 0
        arg_peri = 0
        return CelestialOrbit(self, i, long_asc, arg_peri, eccentricity=0, semimajoraxis=sma)

    def __MakePositionDirty(self):
        self.__last_calculated_time = None
        for c in self.childs:
            c.__MakePositionDirty()

    def setLocalPosition(self, newpos):
        self.__local_pos = newpos

    def getLocalPositionAtTime(self, time=None):
        return self.__local_pos

    def getGlobalPositionAtTime(self, time=None, do_update = False):
        if self.__last_calculated_time == time and time is not None:
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
        except AttributeError as e:
            return position
        if isinstance(obj, CelestialBody):
            par_pos = obj.getGlobalPositionAtTime(time, do_update)
            return position + par_pos
        else:
            return position

    def clearCachedLocation(self):
        self.__last_calculated_time = None

class CelestialStar(CelestialBody):
    def __init__(self, brightness, **kwargs):
        super().__init__(**kwargs)
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
    def __init__(self, SOI_radius, mean_anomaly_at_epoch=0, orbit=None, **kwargs):
        super().__init__(**kwargs)
        if orbit is None:
            self.orbit = None
        else:
            self.orbit = orbit
            p = self.orbit.parent
            p.AppendChild(self)
        self.SOI = SOI_radius
        self.M0 = MeanAnomaly(mean_anomaly_at_epoch)
    def __repr__(self):
        return "CelestialBody_Orbitting(unique_name=%r, mass=%r, mu=%r, radius=%r, sidereal_period=%r, atmosphere_height=%r, low_orbit_height=%r, mean_anomaly_at_epoch = %r, \n\t orbit = %r)" \
               % (self.name, self.mass, self.mu, self.radius, self.sidereal_rotation_period, self.atmosphere_height,
                  self.low_orbit_height, self.M0, self.orbit)

    def getMeanAnomaly(self, time):
        return self.M0 + MeanAnomaly(self.orbit.mean_motion * time)

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

    def getLocalPositionAtTrueAnomaly(self, true_anomaly):
        return self.orbit.TrueAnomalyToPositionVector(true_anomaly)
    def getLocalPositionAtTime(self, time):
        return self.getLocalPositionAtTrueAnomaly(self.getTrueAnomaly(time))

    @property
    def orbit_period(self):
        return self.orbit.OrbitalPeriod()

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
