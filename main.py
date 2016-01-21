
from scipy.optimize import *
from math import *
import matplotlib.pyplot as plt
from celestial_body import *
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
    shadow_time = orbit.TotalMaxTimeInShadow()
    light_time = orbit.OrbitalPeriod() - shadow_time
    energy = shadow_time * use_power
    solar_cell_power = energy / light_time + use_power
    return (energy, solar_cell_power)


def CommSatelliteMinPower(orbit, core, downlink_antenna, constellation_satellite_num, extra_antennas=None,
                          extra_energy=0, distance_safety_margin=0.1):
    extra_antennas_power = 0
    if extra_antennas is not None:
        extra_antennas_power = sum(antenna.power for antenna in extra_antennas)
    min_obj = minimize_scalar(
        lambda x: -orbit.PythagoralDistance(orbit.ToTrueAnomaly(MeanAnomaly(x)),
                                            orbit.ToTrueAnomaly(
                                                MeanAnomaly(x - (2 * pi) / constellation_satellite_num))),
        bounds=(0, pi), method="bounded", options={'maxiter': 1000, 'xatol': np.finfo(float).eps})
    dis = (1 + distance_safety_margin) * -min_obj.fun

    return core.power + extra_energy + extra_antennas_power + downlink_antenna.power  # + communication_antenna.power*constellation_satellite_num

if __name__ == "__main__":
    # d = {t: (0.3, 0.3, 0.3)}
    # for key, value in d.items():
    #    print(key(*value))
    Kerbol = CelestialStar(name="Kerbol", mass=1.7565670E28, mu=1.1723328E18, radius=261600000,
                           sidereal_rotation_period=432000, atmosphere_height=600000, brightness=1)
    kerbin_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=0, longitude_ascending_node=0,
                                  argument_periapsis=0, eccentricity=0, semimajoraxis=13599840256)
    Kerbin = CelestialBody_Orbitting(name="Kerbin", mass=5.2915793e22, mu=3.5316000e15, radius=600000,
                                     sidereal_rotation_period=21599.912, SOI_radius=84159286, atmosphere_height=70000,
                                     mean_anomaly_at_epoch=pi, orbit=kerbin_orbit)

    mun_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=12000000)
    Mun = CelestialBody_Orbitting(name="Mun", mass=9.7600236e20, mu=6.51383980e10,
                                  radius=200000, sidereal_rotation_period=138984.38,
                                  mean_anomaly_at_epoch=1.7, SOI_radius=2429559.1, orbit=mun_orbit)

    minmus_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=radians(6),
                                  longitude_ascending_node=radians(78), argument_periapsis=radians(38), eccentricity=0,
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

    sat_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=Mun.radius +0.0)


    sat = satellite.satellite(sat_orbit)
    sat.addPart(KSP_PART.ControlProbe["probeCoreOcto"])
    for i in range(5):
        sat.addPart(KSP_PART.Antenna["mediumDishAntenna"])
    print(sat.orbit.MaxTimeInShadow(), sat.orbit.TimeInShadow(0) - sat.orbit.MaxTimeInShadow())
    print(sat.getMaxAverageDistance())
    print(sat.getMaxDistance())
    print(sat.power)
    print(satellite.calculate_minimal_electric_storage(sat))
    print(satellite.calculate_minimal_charge_rate(sat))

    #print('{0:.15f}'.format((sat.getMaxDistance() / sat.getMaxAverageDistance())**2))


    print("done")
