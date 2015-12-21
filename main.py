
from scipy.optimize import *
from math import *
import matplotlib.pyplot as plt
from celestial_body import *
from resources import *
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors
import timeit

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
    Kerbol = CelestialBody("Kerbol", mass=1.7565670E28, mu=1.1723328E18, radius=261600000, sidereal_period=432000,
                           atmosphere_height=600000)
    kerbin_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=0, longitude_ascending_node=0,
                                  argument_periapsis=0, eccentricity=0, semimajoraxis=13599840256)
    Kerbin = CelestialBody_Orbitting("Kerbin", mass=5.2915793e22, mu=3.5316000e12, radius=600000,
                                     sidereal_period=21599.912, SOI_radius=84159286, atmosphere_height=70000,
                                     mean_anomaly_at_epoch=pi, orbit=kerbin_orbit)

    mun_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=12000000)
    Mun = CelestialBody_Orbitting("Mun", mass=9.7600236e20, mu=6.51383980e10, radius=200000, sidereal_period=138984.38,
                                  mean_anomaly_at_epoch=1.7, SOI_radius=2429559.1, orbit=mun_orbit)

    minmus_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=radians(6),
                                  longitude_ascending_node=radians(78), argument_periapsis=radians(38), eccentricity=0,
                                  semimajoraxis=47000000)
    Minmus = CelestialBody_Orbitting("Minmus", mass=2.6457897e19, mu=1.7658000e9, radius=60000, sidereal_period=40400,
                                     mean_anomaly_at_epoch=0.9, SOI_radius=2247428.4, orbit=minmus_orbit)

    max_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0, semimajoraxis=Kerbin.SOI)
    max_mun_orbit = CelestialOrbit(celestial_parent_body=Mun, inclination=0, longitude_ascending_node=0,
                                   argument_periapsis=0, eccentricity=0, semimajoraxis=Mun.SOI)
    max_minmus_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=0, longitude_ascending_node=0,
                                      argument_periapsis=0, eccentricity=0, semimajoraxis=Minmus.SOI)

    sat_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=0.9999, semimajoraxis=Minmus.radius + 200000)

    anomalies = np.linspace(0, pi, 200000, endpoint=False)

    fig_time = plt.figure()
    ax_time = fig_time.add_subplot("111")


    theta = np.linspace(0, 2 * pi, 101, endpoint=True)
    v_theta = np.linspace(0,2*pi, 8, endpoint=False)
    sat_number = np.arange(2, 20)
    ecc = np.linspace(0., 0.9, 10)



    e = 0.5
    n = 1000
    a = 1e7
    fig_orbit = plt.figure()
    ax_orbit = fig_orbit.add_subplot(111, projection='3d')
    ax_orbit.set_xlabel("x")
    ax_orbit.set_ylabel("y")
    ax_orbit.set_zlabel("z")
    ax_orbit.set_xlim3d([-1e7,1e7])
    ax_orbit.set_ylim3d([-1e7,1e7])
    ax_orbit.set_zlim3d([-1e7,1e7])
    orbit1 = CelestialOrbit(celestial_parent_body=Kerbin, inclination=pi/4, longitude_ascending_node=0,
                           argument_periapsis=1/4 * pi, eccentricity=e, semimajoraxis=a)
    p1 = orbit1.TrueAnomalyToPositionVector(theta)
    p1_v = orbit1.TrueAnomalyToPositionVector(v_theta)
    v1_v = orbit1.TrueAnomalyToVelocityVector(v_theta)
    orbit2 = CelestialOrbit(celestial_parent_body=Kerbin, inclination=pi/2, longitude_ascending_node=0,
                           argument_periapsis=0, eccentricity=0.5, semimajoraxis=a)
    p2 = orbit2.TrueAnomalyToPositionVector(theta)
    p2_v = orbit2.TrueAnomalyToPositionVector(orbit2.TrueAnomalyFromEccentricAnomaly(v_theta))
    v2_v = orbit2.TrueAnomalyToVelocityVector(orbit2.TrueAnomalyFromEccentricAnomaly(v_theta))
    orbit3 = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0.2, longitude_ascending_node=pi/2,
                           argument_periapsis=pi/6, eccentricity=.3, semimajoraxis=a)
    p3 = orbit3.TrueAnomalyToPositionVector(theta)
    p3_v = orbit3.TrueAnomalyToPositionVector(v_theta)
    v3_v = orbit3.TrueAnomalyToVelocityVector(v_theta)


    print(v2_v)
    i=0
    for _v2 in  v2_v:
        print(v_theta[i], _v2, sqrt(sum(_**2 for _ in _v2)))
        i+=1
    #print(p)
    #ax_orbit.plot(p1[:,0],p1[:,1],p1[:,2])
    ax_orbit.plot(p2[:,0],p2[:,1],p2[:,2])
    ax_orbit.quiver(p2_v[:,0],p2_v[:,1],p2_v[:,2],v2_v[:,0],v2_v[:,1],v2_v[:,2],length=4e6)
    #ax_orbit.plot(p3[:,0],p3[:,1],p3[:,2])

    #ax_time.plot(p[:,0],p[:,1])
    ax_time.plot(p2[:,0],p2[:,1])
    ax_time.plot(p3[:,0],p3[:,1])
    plotSphere(0,0,0,Kerbin.radius,ax_orbit, color='r')


    plt.show()
    #exit()

    sma = np.linspace(Kerbin.radius + Kerbin.atmosphere_height, Kerbin.SOI, 10)
    max_distances = []
    exit()
    for e in ecc:
        orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0,
                               argument_periapsis=0, eccentricity=e, semimajoraxis=a)
        l = lambda x: orbit.PythagoralDistance(
            orbit.ToTrueAnomaly(MeanAnomaly(x)), orbit.ToTrueAnomaly(MeanAnomaly(x - (2 * pi) / n))
        )
        min_obj = minimize_scalar(lambda x: -l(x), bounds=(0, pi), method="bounded",
                                  options={'maxiter': 1000, 'xatol': np.finfo(float).eps})
        print(e, min_obj.x, l(min_obj.x), min_obj)
        max_distances.append(-min_obj.fun / a)
        nu1 = np.array([orbit.ToTrueAnomaly(MeanAnomaly(M)) for M in theta])
        nu2 = np.array([orbit.ToTrueAnomaly(MeanAnomaly(M - (2 * pi) / n)) for M in theta])
        distance = np.array([orbit.PythagoralDistance(n1, n2) for n1, n2 in zip(nu1, nu2)])
        ax_time.plot(nu1 / pi, distance, label="ecc: " + str(e))
        # ax_time.plot(theta,nu)

    ax_time.plot(ecc, max_distances)
    ax_time.set_xlabel("True anomaly (0 = periapsis) $[\pi RAD]$")
    ax_time.set_ylabel("Distance between sats in " +str(n)+" sat constellation $[m]$")
    ax_time.grid()
    ax_time.legend()
    plt.show()



    print("done")
