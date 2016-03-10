from math import pi
import numpy as np
import math as m
from numpy import cos, sin, sqrt, power, square, arctan2, arccos, arcsin, arcsinh, radians, degrees
from scipy.optimize import *
import typing

__author__ = 'Paul'
atan2 = arctan2
acos = arccos
asin = arcsin
asinh = arcsinh



#CelestialBody_Type = typing.TypeVar("CelestialBody_Type", "CB.CelestialBody", np.ndarray)

def _keplerEq(E:float, eccentricity:float, mean_anomaly:float) -> float:
    return E - eccentricity * sin(E) - mean_anomaly


def _keplerEqPrime(E:float, eccentricy:float, mean_anomaly:float) -> float:
    return 1 - eccentricy * cos(E)


def _keplerEqPrime2(E:float, eccentricity:float, mean_anomaly:float) -> float:
    return eccentricity * sin(E)


def _TransformCoordinates(mat:np.ndarray, vec:np.ndarray) -> np.ndarray:
    return mat.dot(vec)


def _VectorRotateOriginAxis(vec:np.ndarray, unit_dir_vector:np.ndarray, theta:float) -> np.ndarray:
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


def _VectorRotateX(vec:np.ndarray, theta:float)->np.ndarray:
    mat = np.array([[1, 0, 0],
                    [0, cos(theta), -sin(theta)],
                    [0, sin(theta), cos(theta)]])
    return mat.dot(vec)


def _VectorRotateY(vec:np.ndarray, theta:float) -> np.ndarray:
    mat = np.array([[cos(theta), 0, -sin(theta)],
                    [0, 1, 0],
                    [sin(theta), 0, cos(theta)]])
    return mat.dot(vec)


def _VectorRotateZ(vec:np.ndarray, theta:float) -> np.ndarray:
    mat = np.array([[cos(theta), -sin(theta), 0],
                    [sin(theta), cos(theta), 0],
                    [0, 0, 1]])
    return mat.dot(vec)


class Anomaly(float):
    """
    Anomaly class
    """

    @classmethod
    def FromDegrees(cls, deg:float) -> "Anomaly":
        """
        :param deg: degrees
        :return: anomaly in radians
        :rtype: Anomaly
        """
        return cls(radians(deg))

    def ToDegrees(self) -> float:
        return degrees(self)

    def ModulateFullCircle(self) -> "Anomaly":
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
    def __init__(self, celestial_parent_body:"CB.CelestialBody", inclination:float,
                 longitude_ascending_node:float, argument_periapsis:float, eccentricity:float=None,
                 semimajoraxis:float=None, semilatusrectum:float=None, periapsis:float=None, apoapsis:float=None):
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
        self._set_orbit_shape(eccentricity, semimajoraxis, semilatusrectum, periapsis, apoapsis) #type:float
        self.parent = celestial_parent_body #type:CB.CelestialBody
        self.__i = inclination #type:float
        self.__long_asc_node = longitude_ascending_node #type:float
        self.__arg_peri = argument_periapsis #type:float
        self.__UpdateConvMats()

    def __repr__(self):
        return "CelestialOrbit(celestial_parent_body=%s, inclination=%r, longitude_ascending_node=%r, argument_periapsis = %r, eccentricity=%r, semimajoraxis = %r)" \
               % (self.parent, self.inclination, self.long_asc_node, self.arg_peri, self.e, self.semi_major_axis)

    def orbital_period(self) -> float:
        """
        T = 2*pi * sqrt(a^3 / mu)
        :return: Period
        :rtype: float
        """
        return 2 * pi * sqrt(self.semi_major_axis ** 3 / self.parent.mu)

    def specific_orbital_energy(self) -> float:
        """
        E = -mu/(2*a)
        :return: Specific Orbital energy [Joule]
        :rtype: float
        """
        return -self.parent.mu / (2 * self.semi_major_axis)

    def vis_viva(self, radius:float) -> float:
        """
        Returns speed given a distance

        :param radius: Radius
        :type radius: float
        :rtype: float
        :return: speed
        """
        return sqrt(2 / radius - 1 / self.semi_major_axis)

    def radius_from_true_anomaly(self, nu:float) -> float:
        """
        :param nu: True anomaly
        :type nu: TrueAnomaly
        :return: Radius
        :rtype: float
        """
        return self.semi_latus_rectum / (1 + self.e * cos(nu))

    def radius_from_eccentric_anomaly(self, E:float) -> float:
        """
        :param E: Eccentric anomaly
        :type E: EccentricAnomaly
        :return: Radius
        :rtype: float
        """
        return self.semi_major_axis * (1 - self.e * cos(E))

    def to_true_anomaly(self, theta:float) -> TrueAnomaly:
        """
        :param theta: Anomaly
        :type theta: Anomaly
        :return: True Anomaly
        :rtype: TrueAnomaly
        """
        if isinstance(theta, EccentricAnomaly):
            return self.true_from_eccentric_anomaly(theta)
        elif isinstance(theta, MeanAnomaly):
            return self.true_from_eccentric_anomaly(self.eccentric_from_mean_anomaly(theta))
        else:
            return TrueAnomaly(theta)

    def true_from_eccentric_anomaly(self, E:float) -> EccentricAnomaly:
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

    def to_eccentric_anomaly(self, theta:float) -> EccentricAnomaly:
        """
        :param theta: Anomaly
        :type theta: Anomaly
        :return: Eccentric Anomaly
        :rtype: EccentricAnomaly
        """
        if isinstance(theta, TrueAnomaly):
            return self.eccentric_from_true_anomaly(theta)
        elif isinstance(theta, MeanAnomaly):
            return self.eccentric_from_mean_anomaly(theta)
        else:
            return EccentricAnomaly(theta)

    def eccentric_from_true_anomaly(self, nu:float) -> EccentricAnomaly:
        """
        :type nu: TrueAnomaly
        :param nu: True anomaly
        :return: Eccentric Anomaly
        :rtype: EccentricAnomaly
        """
        if self.e == 0:
            return nu
        return 2 * atan2(sqrt(1 - self.e) * sin(nu / 2), sqrt(1 + self.e) * cos(nu / 2))

    def _EccentricFromMeanAnomalyNewtonMethod(self, M:float) -> EccentricAnomaly:
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

    def eccentric_from_mean_anomaly(self, M:float) -> EccentricAnomaly:
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

    def to_mean_anomaly(self, theta:float) -> MeanAnomaly:
        """
        :type theta: Anomaly
        :param theta: Anomaly
        :return: Mean Anomaly
        :rtype: MeanAnomaly
        """
        if isinstance(theta, EccentricAnomaly):
            return self.mean_from_eccentric_anomaly(theta)
        elif isinstance(theta, TrueAnomaly):
            return self.mean_from_eccentric_anomaly(self.eccentric_from_true_anomaly(theta))
        else:
            return MeanAnomaly(theta)

    def mean_from_eccentric_anomaly(self, E:float) -> MeanAnomaly:
        """
        :type E: EccentricAnomaly
        :param E: Eccentric anomaly
        :return: Mean Anomaly
        :rtype: MeanAnomaly
        """
        return E - self.e * sin(E)

    def in_shadow(self, solar_direction:float) -> bool:
        """
        :param solar_direction: Solar direction
        :type solar_direction: Anomaly
        :rtype: (TrueAnomaly, TrueAnomaly)
        :return: Start and ending true anomly of the shadow
        """
        avg_theta = (solar_direction + pi) % (2 * pi)
        R_body = self.parent.radius

        def shadow_func(theta):
            r = self.radius_from_true_anomaly(avg_theta + theta)
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

    def time_in_shadow(self, solar_anomaly:float) -> float:
        """
        :param solar_anomaly: Solar direction
        :type solar_anomaly: Anomaly
        :rtype: float
        :return: Start and ending true anomly of the shadow
        """
        theta = self.in_shadow(solar_anomaly)

        if theta[1] >= 2 * pi:
            theta = (theta[0] - 2 * pi, theta[1] - 2 * pi)
        E = tuple(self.eccentric_from_true_anomaly(t) for t in theta)
        M = tuple(self.mean_from_eccentric_anomaly(t) for t in E)

        res = tuple(t / self.mean_motion for t in M)
        return res[1] - res[0]

    def max_time_in_shadow(self) -> float:
        """
        Returns maximum time in shadow.
        Numerical solution, for standard orbits time_in_shadow(0) might be faster
        :return: time
        :rtype: float
        """
        # if self.e == 0:
        return self.time_in_shadow(0)
        # min_obj = minimize_scalar(lambda x: -self.time_in_shadow(x), bounds=(0, pi), method="bounded")
        # return min_obj.x, self.time_in_shadow(min_obj.x)

    def iterative_max_time_in_shadow(self) -> float:
        """
        Returns total max time in shadow
        Iterativelly calls max_time_in_shadow for each orbit until star is found
        :return: time
        :rtype: float
        """
        return sum(i[1] for i in self.generate_shadow_light_time_list())

    def generate_shadow_light_time_list(self) -> float:
        """
        Returns total max time in shadow
        Iterativelly calls max_time_in_shadow for each orbit until star is found
        :return: list of (light, shadow) times
        """
        o = self
        try:
            while not hasattr(o.parent, "brightness"):
                period = o.period
                t = o.max_time_in_shadow()
                yield (period - t, t)
                o = o.parent.orbit
        except AttributeError:
            return
        return

    def _true_anomaly_to_planar_position(self, nu:float) -> np.ndarray:
        r_abs = self.radius_from_true_anomaly(nu)
        return np.array([r_abs * cos(nu), r_abs * sin(nu), 0, 1])

    def pythagoral_distance(self, nu1:float, nu2:float) -> float:
        """
        :type nu1: TrueAnomaly
        :param nu1: Anomaly first object
        :param nu2: Anomaly second object
        :type nu2: TrueAnomaly
        :return: Distance between objects
        :rtype: float
        """
        p1 = self._true_anomaly_to_planar_position(nu1)[:3]
        p2 = self._true_anomaly_to_planar_position(nu2)[:3]
        d = np.linalg.norm(p1 - p2)
        return d

    def true_anomaly_to_position(self, nu:float) -> np.ndarray:
        """
        :param nu: TrueAnomaly
        :type nu: TrueAnomaly
        :return: Vector
        :rtype: np.ndarray
        """
        p_in_plane = self._true_anomaly_to_planar_position(nu)
        p_local = self.__planar_to_local_mat.dot(p_in_plane)
        return p_local[:3]

    def true_anomaly_to_velocity(self, nu:float) -> np.ndarray:
        v_inplane = sqrt(self.parent.mu / self.semi_latus_rectum) * np.array([-sin(nu), self.e + cos(nu), 0, 0])
        v = self.__planar_to_local_mat.dot(v_inplane)
        return v[:3]

    def _set_orbit_shape(self, eccentricity:float=None, semimajoraxis:float=None, semilatusrectum:float=None,
                         periapsis:float=None, apoapsis:float=None):
        t = [eccentricity is not None, semimajoraxis is not None, semilatusrectum is not None, periapsis is not None,
             apoapsis is not None]
        str = "Incorrect amount of orbital variables"
        if sum(t) != 2:
            raise TypeError(str)
        if eccentricity is not None:
            self.e = eccentricity
            if semilatusrectum is not None:
                self.semi_latus_rectum = semilatusrectum
            elif semimajoraxis is not None:
                self.semi_major_axis = semimajoraxis
            elif apoapsis is not None:
                self.semi_major_axis = semimajoraxis
            elif periapsis is not None:
                self.semi_latus_rectum = periapsis * (eccentricity + 1)
            else:
                raise TypeError(str)
        elif semilatusrectum is not None:
            if periapsis is not None:
                self.e = semilatusrectum / periapsis - 1
                raise TypeError(str)
        elif semimajoraxis is not None:
            self.semi_major_axis = semimajoraxis
            if apoapsis is not None:
                periapsis = 2 * semimajoraxis - apoapsis
            if periapsis is not None:
                self.e = 1 - periapsis / semimajoraxis
            else:
                raise TypeError(str)

    def create_tree_branch(self, ancestor_body:"CB.CelestialBody"=None) -> typing.Generator["CelestialOrbit",None,None]:
        orbit = self
        yield self
        while ancestor_body is None or orbit.parent is not ancestor_body:
            try:
                orbit = orbit.parent.orbit
                yield orbit
            except AttributeError:
                break

    def get_normal_to_orbital_plane(self) -> np.ndarray:
        # vec1 = np.array([cos(self.long_asc_node), sin(self.long_asc_node),0])
        # vec2 = np.array([-sin(self.long_asc_node)*cos(self.inclination), cos(self.long_asc_node)*cos(self.inclination), sin(self.inclination)])
        # v = np.cross(vec1, vec2)
        sin_incl = sin(self.inclination)
        v = np.array([sin_incl * sin(self.long_asc_node), -sin_incl * cos(self.long_asc_node), cos(self.inclination)])
        return v

    def get_orbital_plane_constants(self) -> tuple:
        n = self.get_normal_to_orbital_plane()
        p = self.parent.get_global_position_at_time()
        t = -n * p
        return tuple(n) + (np.sum(t),)

    def project_local_point_on_orbital_plane(self, point: np.ndarray) -> np.ndarray:
        return self._project_local_point_on_orbital_plane(point, self.get_normal_to_orbital_plane())
    def _project_local_point_on_orbital_plane(self, point: np.ndarray, normal:np.ndarray) -> np.ndarray:
        return point - np.dot(point, normal) * normal

    def to_planar_coordinates(self, local_3d_on_plane: np.ndarray) -> np.ndarray:
        return self.__local_to_planar_mat[:3, :3].dot(local_3d_on_plane)[:2]

    def to_local_coordinates(self, point_2d_in_plane:np.ndarray) -> np.ndarray:
        p = np.array([point_2d_in_plane[0], point_2d_in_plane[1], 0, 1])
        return self.__planar_to_local_mat.dot(p)[:3]

    def to_global_position_at_time(self, local_coordinates, time=None):
        return local_coordinates + self.parent.get_global_position_at_time(time)

    def get_distance(self, root_point: np.ndarray, true_anomaly: TrueAnomaly) -> float:
        p = self.true_anomaly_to_position(true_anomaly)
        parent_pos = self.parent.get_global_position_at_time()
        d = np.linalg.norm(p + parent_pos - root_point)
        return d

    def getMaxDistance_outdated(self, ancestor_body:"CB.CelestialBody") -> float:
        if self.parent == ancestor_body:
            return self.apoapsis_distance
        # l = list(self.create_tree_branch(ancestory_body))
        root_point = ancestor_body.get_global_position_at_time()

        def fun_negative_distance(x, root_point):
            return -self.get_distance(root_point, x)

        guesses = (0, pi / 2, pi, 1.5 * pi)
        args = (root_point,)
        guess = max((self.get_distance(root_point, guess), guess) for guess in guesses)[1]
        bracket = (guess - pi / 2, guess, guess + pi / 2)
        res = minimize_scalar(fun_negative_distance,
                              bracket=bracket, method="brent", args=args)
        # print(res.x, fun_negative_distance(res.x))
        return -res.fun

    def get_total_max_distance(self, ancestor_body:"typing.Union[CB.CelestialBody,np.ndarray]",
                               eps:float=3*np.finfo(float).eps) -> float:
        try:
            pos = ancestor_body.get_global_position_at_time()
        except AttributeError:
            pos = ancestor_body
            orbit_list = list(self.create_tree_branch(None))
        else:
            if self.parent is ancestor_body:
                return self.apoapsis_distance
            orbit_list = list(self.create_tree_branch(ancestor_body))
        orbit_list.reverse()
        return orbit_list[0]._get_total_max_distance(pos, pos, orbit_list[1:], eps)

    def _get_total_max_distance(self, root_point:np.ndarray,
                                parent_position:np.ndarray,
                                orbit_list:typing.List["CelestialOrbit"],
                                eps: float = 3*np.finfo(float).eps) -> float:
        other = orbit_list[0]

        def calculate_distance_(true_anomaly):
            p = parent_position + self.true_anomaly_to_position(true_anomaly)
            t = p - root_point
            base_distance = m.sqrt(t.dot(t))
            t = base_distance * eps
            sma = other.semi_major_axis
            if sma <= t:
                return base_distance
            # recursive optimization. (for each position find the optimal solution of the child's orbit)
            if len(orbit_list) == 1:
                return other._get_max_distance_non_recursive(root_point, p, eps)
            else:
                return other._get_total_max_distance(root_point, p, orbit_list[1:], eps)

        res = minimize_scalar(lambda x: -calculate_distance_(x), bounds=(0, 2 * pi - eps), method="Bounded", tol=eps)
        return -res.fun

    def get_max_distance_non_recursive(self, ancestor_body:"CB.CelestialBody",
                                       eps:float=3*np.finfo(float).eps) -> float:
        if self.parent == ancestor_body:
            return self.apoapsis_distance
        root_point = ancestor_body.get_global_position_at_time()
        parent_position = self.parent.get_global_position_at_time()
        return self._get_max_distance_non_recursive(root_point, parent_position, eps)

    def _get_max_distance_non_recursive(self, root_point:np.ndarray, parent_position:np.ndarray,
                                        eps:float=3*np.finfo(float).eps) -> float:
        local_root = root_point - parent_position
        plane_normal = self.get_normal_to_orbital_plane()
        distance_to_plane = plane_normal.dot(local_root)
        proj_root = self._project_local_point_on_orbital_plane(local_root, plane_normal)
        if self.e == 0:
            in_plane_distance_squared = (m.sqrt(proj_root.dot(proj_root)) + self.semi_major_axis)**2
        else:
            ellipse_root = self.to_planar_coordinates(proj_root)
            # origin transform from focal point to center
            ellipse_root[0] += self.e * self.semi_major_axis
            if abs(ellipse_root[1] / ellipse_root[0]) <= eps:
                # on long axis, or center ellipse, quick way to calculate
                in_plane_distance_squared = (abs(ellipse_root[0]) + self.semi_major_axis)**2
            else:
                if ellipse_root[0] >= 0:
                    if ellipse_root[1] > 0:
                        left_bound = pi
                        right_bound = 1.5 * pi - eps
                    else:
                        left_bound = 0.5 * pi + eps
                        right_bound = pi
                else:
                    if ellipse_root[1] > 0:
                        left_bound = 1.5 * pi + eps
                        right_bound = 2 * pi
                    else:
                        left_bound = 0
                        right_bound = 0.5 * pi - eps

                a = self.semi_major_axis
                b = self.semi_minor_axis

                def fun_planar_distance(x, a, b, p):
                    return sin(2.0 * x) / 2 * (b ** 2 - a ** 2) + \
                           p[0] * a * sin(x) - p[1] * b * cos(x)

                args = (a, b, ellipse_root)
                try:
                    theta = brentq(fun_planar_distance, left_bound, right_bound,
                                   args=args, rtol=eps * 2)
                except ValueError:
                    v_left = fun_planar_distance(left_bound, *args)
                    v_right = fun_planar_distance(right_bound, *args)
                    if v_left == v_right:
                        theta = (left_bound + right_bound) / 2
                    elif v_left < v_right:
                        theta = left_bound
                    else:
                        theta = right_bound
                in_plane_distance_squared = (a * cos(theta) - ellipse_root[0]) ** 2 +\
                                            (b * sin(theta) - ellipse_root[1]) ** 2

        return m.sqrt(distance_to_plane ** 2 + in_plane_distance_squared)

    def get_average_distance_current_time(self) -> float:
        body = self.parent
        if hasattr(body, "brightness"):
            return self.semi_major_axis
        while not hasattr(body.orbit.parent, "brightness"):
            body = body.orbit.parent
        return body.get_distance_central_body()

    @property
    def inclination(self) -> float:
        return self.__i

    @inclination.setter
    def inclination(self, value:float):
        self.__i = value
        self.__UpdateConvMats()

    @property
    def long_asc_node(self) -> float:
        return self.__long_asc_node

    @long_asc_node.setter
    def long_asc_node(self, v:float):
        self.__long_asc_node = v
        self.__UpdateConvMats()

    @property
    def arg_peri(self) -> float:
        return self.__arg_peri

    @arg_peri.setter
    def arg_peri(self, v:float):
        self.__arg_peri = v
        self.__UpdateConvMats()

    @property
    def semi_major_axis(self) -> float:
        return self.semi_latus_rectum / (1 - self.e ** 2)

    @semi_major_axis.setter
    def semi_major_axis(self, value:float):
        self.semi_latus_rectum = value * (1 - self.e ** 2)

    @property
    def semi_minor_axis(self) -> float:
        return self.semi_latus_rectum / sqrt(1 - self.e ** 2)

    @semi_minor_axis.setter
    def semi_minor_axis(self, value:float):
        self.semi_latus_rectum = value * sqrt(1 - self.e ** 2)

    @property
    def mean_motion(self) -> float:
        return sqrt(self.parent.mu / self.semi_major_axis ** 3.)

    @property
    def period(self) -> float:
        return 2 * pi * sqrt(self.semi_major_axis ** 3 / self.parent.mu)

    @property
    def periapsis_distance(self) -> float:
        return self.semi_latus_rectum / (1 + self.e)

    @property
    def apoapsis_distance(self) -> float:
        return self.semi_latus_rectum / (1 - self.e)

    def __UpdateConvMats(self):
        t = np.array([[cos(self.long_asc_node), -sin(self.long_asc_node), 0, 0],
                      [sin(self.long_asc_node), cos(self.long_asc_node), 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
        t = t.dot([[1, 0, 0, 0],
                   [0, cos(self.inclination), -sin(self.inclination), 0],
                   [0, sin(self.inclination), cos(self.inclination), 0],
                   [0, 0, 0, 1]])
        t = t.dot([[cos(self.arg_peri), -sin(self.arg_peri), 0, 0],
                   [sin(self.arg_peri), cos(self.arg_peri), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])
        self.__planar_to_local_mat = t
        self.__local_to_planar_mat = t.T
