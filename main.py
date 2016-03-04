import random
from scipy.optimize import *
from math import *
import matplotlib.pyplot as plt
from celestial_body import *
from celestial_orbit import *
from ksp_resources import *
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors
import timeit
import scipy as sp
import satellite
import ksp_part_resource as KSP_PART
import ksp_part_remotetech
import ksp_part_squad
import copy
import itertools

def FuelFraction(dV, I_sp, g0 = 9.81):
    # dv = I_sp*g0 * ln(M0/M_dry)
    # M0/M_dry = exp(dV/v_e)
    # M0/(M0 - M_f) = exp(dV/v_e)
    # (M0 - M_f)/M0 = 1/exp(dV/v_e)
    # 1 - FF = exp(-dV/v_e)
    # FF = 1 - exp(-dV/v_e)
    return 1 - math.exp(- dV /(I_sp * g0))

def FuelMass(dV, payload_mass, I_sp, g0 = 9.81):
    # ff = M_fuel / M0 = M_fuel / (M_dry + M_fuel)
    # M_dry = M_fuel* (1 - ff)
    return payload_mass/ (1 - FuelFraction(dV, I_sp, g0))


def ForceTotal(I_sp, m_dot, g_0=9.81):
    return I_sp * m_dot * g_0


def MassFlow(F_thrust, I_sp, g_0=9.81):
    return F_thrust / (I_sp * g_0)


def TotalFuelMass(MassFlow, BurnTime):
    return MassFlow * BurnTime


def BurnTime(dV, mass_total, F, I_sp, g0=9.81):
    m_dot = MassFlow(F, I_sp, g0)
    FF = FuelFraction(dV, I_sp, g0)
    mass_fuel = mass_total * FF

    return (mass_fuel / m_dot, mass_fuel, m_dot, FF)


def delta_v(payload_mass, fuel_mass, I_sp, g0=9.81):
    return I_sp * g0 * math.log((fuel_mass + payload_mass) / payload_mass)


def Satellite_ElectricalDesign(orbit, use_power):
    shadow_time = orbit.total_max_time_in_shadow()
    light_time = orbit.orbital_period() - shadow_time
    energy = shadow_time * use_power
    solar_cell_power = energy / light_time + use_power
    return (energy, solar_cell_power)


def CommSatelliteMinPower(orbit, core, downlink_antenna, constellation_satellite_num, extra_antennas=None,
                          extra_energy=0, distance_safety_margin=0.1):
    extra_antennas_power = 0
    if extra_antennas is not None:
        extra_antennas_power = sum(antenna.power for antenna in extra_antennas)
    min_obj = minimize_scalar(
        lambda x: -orbit.pythagoral_distance(orbit.to_true_anomaly(MeanAnomaly(x)),
                                            orbit.to_true_anomaly(
                                                MeanAnomaly(x - (2 * pi) / constellation_satellite_num))),
        bounds=(0, pi), method="bounded", options={'maxiter': 1000, 'xatol': np.finfo(float).eps})
    dis = (1 + distance_safety_margin) * -min_obj.fun

    return core.power + extra_energy + extra_antennas_power + downlink_antenna.power  # + communication_antenna.power*constellation_satellite_num



class TestObj(object):
    def __init__(self, payload=None):
        self.__all_parts = []
        self.payload = payload
        self.structural_mass = 0

    @property
    def empty_mass(self):
        return self.structural_mass + sum(part.mass for part in self.__all_parts if hasattr(part, "mass"))

    @property
    def mass(self):
        mass = self.empty_mass
        try:
            mass += self.payload.mass
        except AttributeError:
            pass
        return mass

if __name__ == "__main__":
    Kerbol = CelestialStar(name="Kerbol", mass=1.7565670E28, mu=1.1723328E18, radius=261600000,
                           sidereal_rotation_period=432000, atmosphere_height=600000, brightness=1)
    a = 0.1*pi
    kerbin_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=0, longitude_ascending_node=a,
                                  argument_periapsis=0, eccentricity=0.000, semimajoraxis=13599840256)
    Kerbin = CelestialBody_Orbitting(name="Kerbin", mass=5.2915793e22, mu=3.5316000e15, radius=600000,
                                     sidereal_rotation_period=21599.912, SOI_radius=84159286, atmosphere_height=70000,
                                     mean_anomaly_at_epoch=0.5, orbit=kerbin_orbit)

    mun_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=12000000)
    Mun = CelestialBody_Orbitting(name="Mun", mass=9.7600236e20, mu=6.51383980e10,
                                  radius=200000, sidereal_rotation_period=138984.38,
                                  mean_anomaly_at_epoch=1.7, SOI_radius=2429559.1, orbit=mun_orbit)

    minmus_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=radians(6),
                                  longitude_ascending_node=radians(78), argument_periapsis=radians(38), eccentricity=0.0,
                                  semimajoraxis=47000000)
    Minmus = CelestialBody_Orbitting(name="Minmus", mass=2.6457897e19, mu=1.7658000e9,
                                     radius=60000, sidereal_rotation_period=40400,
                                     mean_anomaly_at_epoch=0.9, SOI_radius=2247428.4, orbit=minmus_orbit)

    max_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=Kerbin.SOI)
    max_mun_orbit = CelestialOrbit(celestial_parent_body=Mun, inclination=0, longitude_ascending_node=0,
                                   argument_periapsis=0, eccentricity=0, semimajoraxis=Mun.SOI)
    max_minmus_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=0, longitude_ascending_node=0,
                                      argument_periapsis=0, eccentricity=0, semimajoraxis=Minmus.SOI)

    sat_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=pi, longitude_ascending_node=pi/3,
                               argument_periapsis=pi, eccentricity=0.99, semimajoraxis=Minmus.SOI)


    techs = [
        "basicScience",
        "start",
        "electrics",
        "advElectrics",
        "largeElectrics",
        "experimentalElectrics",
        "specializedElectrics"
    ]

    #print(Kerbin.orbit_period/2, sat_orbit.period, (Kerbin.orbit_period/2)/sat_orbit.period)
    #x = np.linspace(0, kerbin_orbit.period/2, 100)
    #print(x[1] - x[0])
    y = []
    y2=[]
    y3 = []
    #fig = plt.figure()
    #ax = fig.add_subplot("111")
    #for i in np.nditer(x):
    #    Kerbol.updateTimeWholeSystem(i)
    #    y.append(minmus_orbit.get_max_distance_non_recursive(Kerbol))
    #    y2.append(sat_orbit.get_max_distance_non_recursive(Kerbol))
    #    y3.append(minmus_orbit.get_average_distance_current_time())

    #ax.plot(x,np.array(y3))
    #ax.plot(x,np.array(y))
    #ax.plot(x,np.array(y2))
    print("----------")
    for i in range(1):
        sat_orbit.get_total_max_distance(Kerbol, eps=0.00001)
    #ax.set_ylim([0, ax.get_ylim()[1]])
    #plt.show()
    exit()




    sat = satellite.satellite(orbit=sat_orbit)

    sat.addPart(KSP_PART.ControlProbe["probeCoreOcto"])
    for i in range(5):
        sat.addPart(KSP_PART.Antenna["mediumDishAntenna"])
    print(sat.orbit.max_time_in_shadow(), sat.orbit.time_in_shadow(0) - sat.orbit.max_time_in_shadow())
    print(sat.getMaxAverageDistance())
    print(sat.getMaxDistance())
    #KSP_PART.BatteryPack["batteryPack"].price = 90
    #KSP_PART.BatteryPack["batteryBank"].mass *= 1000

    def check_battery(l, r, useful_storage = np.inf):
        left_power = min(l.electric_charge, useful_storage)
        right_power = min(r.electric_charge, useful_storage)
        left_v = (left_power/l.mass, left_power/l.price, left_power)
        right_v = (right_power/r.mass, right_power/r.price, right_power)
        t = any(ll > rr for ll,rr in zip(left_v, right_v))
        return t

    def check_solar_cell(l, r, distance, useful_power = np.inf):
        left_power = min(l.getChargeRate(distance), useful_power)
        right_power = min(r.getChargeRate(distance), useful_power)
        left_v = (left_power/l.mass, left_power/l.price, left_power)
        right_v = (right_power/r.mass, right_power/r.price, right_power)
        t = any(ll > rr for ll,rr in zip(left_v, right_v))
        return t


    power_requirements = satellite.calculate_minimal_charge_rate(sat)
    print(power_requirements)
    all_copies = satellite.generate_all_solar_cells(sat, techs, power_requirements, first_order_check=check_solar_cell)

    #for i in range(1):
    #    energy_storage_requirements = satellite.calculate_minimal_electric_storage(sat)
    #    #print(energy_storage_requirements)
    #    all_copies = satellite.generate_all_battery(sat,techs, energy_storage_requirements, check_battery)
    print(len(all_copies))
    for s in all_copies:
        print({"mass:": s.mass, "price":s.price})
        for p in s._all_parts:
            print(" - ", p.title)

    print("done")
