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
    return 1 - math.exp(-dV / (I_sp * g0))

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
    shadow_time = orbit.iterative_max_time_in_shadow()
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

    moho_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=radians(7), longitude_ascending_node=radians(70),
                                  argument_periapsis=radians(15), eccentricity=0.2, semimajoraxis=5263138304)
    Moho = CelestialBody_Orbitting(name="Moho", mass=2.5263617e21, mu=1.6860938e11, radius=250000,
                                     sidereal_rotation_period=1210000, SOI_radius=9646663, atmosphere_height=0,
                                     mean_anomaly_at_epoch=3.14, orbit=moho_orbit)
    eve_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=radians(2.1), longitude_ascending_node=radians(15),
                                  argument_periapsis=radians(0), eccentricity=0.01, semimajoraxis=9832684544)
    Eve = CelestialBody_Orbitting(name="Eve", mass=1.2244127e23, mu=8.1717302e12, radius=700000,
                                     sidereal_rotation_period=81661.857, SOI_radius=85109365, atmosphere_height=90000,
                                     mean_anomaly_at_epoch=3.14, orbit=eve_orbit)
    gilly_orbit = CelestialOrbit(celestial_parent_body=Eve, inclination=radians(12), longitude_ascending_node=radians(80),
                                  argument_periapsis=radians(10), eccentricity=0.55, semimajoraxis=31500000)
    Gilly = CelestialBody_Orbitting(name="Gilly", mass=1.2420512e17, mu=8289449.8, radius=13000,
                                     sidereal_rotation_period=28255, SOI_radius=126123.27, atmosphere_height=0,
                                     mean_anomaly_at_epoch=0.9, orbit=gilly_orbit)

    kerbin_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=0, longitude_ascending_node=0,
                                  argument_periapsis=0, eccentricity=0.000, semimajoraxis=13599840256)
    Kerbin = CelestialBody_Orbitting(name="Kerbin", mass=5.2915793e22, mu=3.5316000e15, radius=600000,
                                     sidereal_rotation_period=21549.425, SOI_radius=84159286, atmosphere_height=70000,
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
    duna_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=radians(0.06), longitude_ascending_node=radians(135.5),
                                  argument_periapsis=radians(0), eccentricity=0.05, semimajoraxis=20726155264)
    Duna = CelestialBody_Orbitting(name="Duna", mass=4.5154812e12, mu=3.0136321e11, radius=320000,
                                     sidereal_rotation_period=65517.859, SOI_radius=47921949, atmosphere_height=50000,
                                     mean_anomaly_at_epoch=3.14, orbit=duna_orbit)
    ike_orbit = CelestialOrbit(celestial_parent_body=Duna, inclination=radians(0.2), longitude_ascending_node=radians(0),
                                  argument_periapsis=radians(0), eccentricity=0.03, semimajoraxis=3200000)
    Ike = CelestialBody_Orbitting(name="Ike", mass=2.7821949e20, mu=1.8568369e10, radius=130000,
                                     sidereal_rotation_period=65517.862, SOI_radius=1049598.9, atmosphere_height=0,
                                     mean_anomaly_at_epoch=1.7, orbit=ike_orbit)

    dres_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=radians(5), longitude_ascending_node=radians(280),
                                  argument_periapsis=radians(90), eccentricity=0.14, semimajoraxis=40839348203)
    Dres = CelestialBody_Orbitting(name="Dres", mass=3.2191322e20, mu=2.1484489e10, radius=138000,
                                     sidereal_rotation_period=34800, SOI_radius=32832840, atmosphere_height=0,
                                     mean_anomaly_at_epoch=3.14, orbit=dres_orbit)
    jool_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=radians(1.304), longitude_ascending_node=radians(52),
                                  argument_periapsis=radians(0), eccentricity=0.05, semimajoraxis=68773560320)
    Jool = CelestialBody_Orbitting(name="Jool", mass=4.2332635e24, mu=2.8252800e14, radius=6000000,
                                     sidereal_rotation_period=36000, SOI_radius=2.4559852e9, atmosphere_height=200000,
                                     mean_anomaly_at_epoch=0.1, orbit=jool_orbit)
    laythe_orbit = CelestialOrbit(celestial_parent_body=Jool, inclination=radians(0), longitude_ascending_node=radians(0),
                                  argument_periapsis=radians(0), eccentricity=0, semimajoraxis=27184000)
    Laythe = CelestialBody_Orbitting(name="Laythe", mass=2.9397663e22, mu=1.9620000e12, radius=500000,
                                     sidereal_rotation_period=52980.879, SOI_radius=3723645.8, atmosphere_height=50000,
                                     mean_anomaly_at_epoch=3.14, orbit=laythe_orbit)
    vall_orbit = CelestialOrbit(celestial_parent_body=Jool, inclination=radians(0), longitude_ascending_node=radians(0),
                                  argument_periapsis=radians(0), eccentricity=0, semimajoraxis=43152000)
    Vall = CelestialBody_Orbitting(name="Vall", mass=3.1088028e21, mu=2.0748150e11, radius=300000,
                                     sidereal_rotation_period=105962.09, SOI_radius=2406401.4, atmosphere_height=0,
                                     mean_anomaly_at_epoch=0.9, orbit=vall_orbit)
    tylo_orbit = CelestialOrbit(celestial_parent_body=Jool, inclination=radians(0.025), longitude_ascending_node=radians(0),
                                  argument_periapsis=radians(0), eccentricity=0, semimajoraxis=68500000)
    Tylo = CelestialBody_Orbitting(name="Tylo", mass=4.2332635e22, mu=2.852800e12, radius=600000,
                                     sidereal_rotation_period=211926.36, SOI_radius=10856518, atmosphere_height=0,
                                     mean_anomaly_at_epoch=3.14, orbit=tylo_orbit)
    bop_orbit = CelestialOrbit(celestial_parent_body=Jool, inclination=radians(15), longitude_ascending_node=radians(10),
                                  argument_periapsis=radians(25), eccentricity=0.24, semimajoraxis=128500000)
    Bop = CelestialBody_Orbitting(name="Bop", mass=3.7261536e19, mu=2.4868349e9, radius=65000,
                                     sidereal_rotation_period=544507.43, SOI_radius=1221060.9, atmosphere_height=0,
                                     mean_anomaly_at_epoch=0.9, orbit=bop_orbit)
    pol_orbit = CelestialOrbit(celestial_parent_body=Jool, inclination=radians(4.25), longitude_ascending_node=radians(2),
                                  argument_periapsis=radians(15), eccentricity=0.17, semimajoraxis=179890000)
    Pol = CelestialBody_Orbitting(name="Pol", mass=1.0813636e19, mu=7.2170208e8, radius=44000,
                                     sidereal_rotation_period=901902.62, SOI_radius=1042138.9, atmosphere_height=0,
                                     mean_anomaly_at_epoch=0.9, orbit=pol_orbit)

    eeloo_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=radians(6.15), longitude_ascending_node=radians(50),
                                  argument_periapsis=radians(260), eccentricity=0.26, semimajoraxis=90118820000)
    Eeloo = CelestialBody_Orbitting(name="Eeloo", mass=1.1149358e21, mu=7.4410815e10, radius=210000,
                                     sidereal_rotation_period=19460, SOI_radius=1.1908294e8, atmosphere_height=0,
                                     mean_anomaly_at_epoch=3.14, orbit=eeloo_orbit)

    max_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=Kerbin.SOI)
    max_mun_orbit = CelestialOrbit(celestial_parent_body=Mun, inclination=0, longitude_ascending_node=0,
                                   argument_periapsis=0, eccentricity=0, semimajoraxis=Mun.SOI)
    max_minmus_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=0, longitude_ascending_node=0,
                                      argument_periapsis=0, eccentricity=0, semimajoraxis=Minmus.SOI)

    sat_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=pi, longitude_ascending_node=pi/3,
                               argument_periapsis=pi, eccentricity=0.0, semimajoraxis=Minmus.SOI)


    techs = [
        "basicScience",
        "start",
        "electrics",
        "advElectrics",
        "largeElectrics",
        "experimentalElectrics",
        "specializedElectrics"
    ]

    sat = satellite.Satellite(orbit=sat_orbit, name="sat1") #type:satellite.Satellite
    print("here")
    sat.addPart(KSP_PART.ControlProbe["probeCoreOcto"])
    for t_sat in range(5):
        sat.addPart(KSP_PART.Antenna["mediumDishAntenna"])
    print(sat.orbit.max_time_in_shadow(), sat.orbit.time_in_shadow(0) - sat.orbit.max_time_in_shadow())
    print(sat.average_max_distance_to_star())
    print(sat.max_distance_to_star())
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
        print(l.title, r.title)
        print(left_power, right_power)
        print(left_v, right_v)
        print([ll > rr for ll,rr in zip(left_v, right_v)])
        t = any(ll > rr for ll,rr in zip(left_v, right_v))
        return t


    power_requirements = satellite.calculate_minimal_charge_rate(sat)
    print(power_requirements)
    energy_storage_requirements = sat.calculate_minimal_electric_storage()
    energy_storage_requirements = 3900
    r = satellite.generate_all_battery(techs, energy_storage_requirements, check_battery)
    for t_sat in r:
        print("----------")
        m = sum(i.mass for i in t_sat)
        p = sum(i.price for i in t_sat)
        print("price:", p, "Mass", m)
        for j in t_sat:
            print(j.name)

    print("---------------------------")
    all_copies = satellite.generate_all_solar_cells(techs, power_requirements,
                                                    distance=sat.max_distance_to_star(),
                                                    first_order_check=check_solar_cell)
    for t_sat in all_copies:
        print("----------")
        m = sum(i.mass for i in t_sat)
        p = sum(i.price for i in t_sat)
        print("price:", p, "Mass", m)
        for j in t_sat:
            print(j.title)
    #for i in range(1):
    #    energy_storage_requirements = satellite.calculate_minimal_electric_storage(sat)
    #    #print(energy_storage_requirements)
    #    all_copies = satellite.generate_all_battery(sat,techs, energy_storage_requirements, check_battery)


    print("done")
