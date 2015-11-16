

from scipy.optimize import *
from math import *
import matplotlib.pyplot as plt
from celestial_body import *
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
import matplotlib.colors

def FuelFraction(dV, I_sp, g0 = 9.81):
    #dv = I_sp*g0 * ln(M0/M_dry)
    #M0/M_dry = exp(dV/v_e)
    #M0/(M0 - M_f) = exp(dV/v_e)
    #(M0 - M_f)/M0 = 1/exp(dV/v_e)
    #1 - FF = exp(-dV/v_e)
    #FF = 1 - exp(-dV/v_e)
    return 1 - math.exp(-dV/(I_sp * g0))

def FuelMass(dV, payload_mass, I_sp, g0 = 9.81):
    #ff = M_fuel / M0 = M_fuel / (M_dry + M_fuel)
    #M_dry = M_fuel* (1 - ff)
    return payload_mass/(1 - FuelFraction(dV, I_sp, g0))

def ForceTotal(I_sp, m_dot, g_0 = 9.81):
    return I_sp * m_dot * g_0

def MassFlow(F_thrust, I_sp, g_0 = 9.81):
    return F_thrust/(I_sp * g_0)

def TotalFuelMass(MassFlow, BurnTime):
    return MassFlow * BurnTime



def BurnTime(dV, mass_total, F, I_sp, g0 = 9.81):
    m_dot = MassFlow(F, I_sp, g0)
    FF = FuelFraction(dV, I_sp, g0)
    mass_fuel = mass_total * FF

    return (mass_fuel / m_dot, mass_fuel, m_dot, FF)





def delta_v(payload_mass, fuel_mass,I_sp, g0 = 9.81 ):
    return I_sp * g0 * math.log((fuel_mass + payload_mass) / payload_mass)





def time_calc(theta_avg, O):
    theta = O.InShadow(theta_avg)
    return theta[1] - theta[0]
    if theta[1] >= 2*pi:
        theta = (theta[0] - 2*pi, theta[1] - 2*pi)


    E = tuple(O.EccentricFromTrueAnomaly(t) for t in theta)
    M = tuple(O.MeanFromEccentricAnomaly(t) for t in E)

    t = tuple(t/O.mean_motion for t in M)
    #print(theta,E,M, t)
    return t[1]-t[0]




if __name__ == "__main__":
    Kerbol = CelestialBody("Kerbol", mass=1.7565670E28, mu = 1.1723328E18, radius=261600000, sidereal_period=432000, atmosphere_height=600000)
    kerbin_orbit = CelestialOrbit(celestial_parent_body=Kerbol, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = 0, semimajoraxis = 13599840256)
    Kerbin = CelestialBody_Orbitting("Kerbin", mass=5.2915793e22, mu=3.5316000e12, radius=600000, sidereal_period=21599.912, SOI_radius=84159286, atmosphere_height=70000, mean_anomaly_at_epoch=pi, orbit=kerbin_orbit)

    mun_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = 0, semimajoraxis = 12000000)
    Mun = CelestialBody_Orbitting("Mun", mass=9.7600236e20, mu=6.51383980e10, radius=200000, sidereal_period=138984.38, mean_anomaly_at_epoch=1.7,SOI_radius=2429559.1,  orbit=mun_orbit)

    minmus_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=radians(6), longitude_ascending_node=radians(78), argument_periapsis=radians(38), eccentricity = 0, semimajoraxis = 47000000)
    Minmus = CelestialBody_Orbitting("Minmus", mass=2.6457897e19, mu=1.7658000e9, radius=60000, sidereal_period=40400, mean_anomaly_at_epoch=0.9, SOI_radius=2247428.4, orbit=minmus_orbit)

    max_orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = 0, semimajoraxis = Kerbin.SOI)
    max_mun_orbit = CelestialOrbit(celestial_parent_body=Mun, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = 0, semimajoraxis = Mun.SOI)
    max_minmus_orbit = CelestialOrbit(celestial_parent_body=Minmus, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = 0, semimajoraxis = Minmus.SOI)

    #print("maximum time in shadow: ", max_orbit.TotalMaxTimeInShadow()/3600)
    #print("maximum time in shadow: ", max_mun_orbit.TotalMaxTimeInShadow()/3600)
    #print("maximum time in shadow: ", max_minmus_orbit.TotalMaxTimeInShadow()/3600)



    ecc = np.linspace(0,0.8,9)
    ecc = np.append(ecc,0.85)
    ecc = np.append(ecc,np.linspace(0.9,0.94,5))




    # O = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = 0.94, semimajoraxis = 1e7)
    # theta = O.InShadow(pi)
    # if theta[1] >= 2*pi:
    #     theta = (theta[0] - 2*pi, theta[1] - 2*pi)
    # E = tuple(O.EccentricFromTrueAnomaly(t) for t in theta)
    # M = tuple(O.MeanFromEccentricAnomaly(t) for t in E)
    # print(theta, E, (E[1]-E[0])/pi, M, (M[1]-M[0])/pi)
    #
    # theta = O.InShadow(pi/2)
    # if theta[1] >= 2*pi:
    #     theta = (theta[0] - 2*pi, theta[1] - 2*pi)
    # E = tuple(O.EccentricFromTrueAnomaly(t) for t in theta)
    # M = tuple(O.MeanFromEccentricAnomaly(t) for t in E)
    # print(theta, E, (E[1]-E[0])/pi, M, (M[1]-M[0])/pi)
    #
    # theta = O.InShadow(0)
    # if theta[1] >= 2*pi:
    #     theta = (theta[0] - 2*pi, theta[1] - 2*pi)
    # E = tuple(O.EccentricFromTrueAnomaly(t) for t in theta)
    # M = tuple(O.MeanFromEccentricAnomaly(t) for t in E)
    # print(theta, E, (E[1]-E[0])/pi, M, (M[1]-M[0])/pi)
    #
    # for e in ecc:
    #     O = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = e, semimajoraxis = 1e7)
    #     E_x = np.linspace(0,2*pi)
    #     M_y = np.array([O.MeanFromEccentricAnomaly(E)for E in E_x])
    #     plt.plot(E_x, M_y, label="e"+str(e))
    # plt.show()
    #
    #
    # exit()


    fig_time = plt.figure()
    ax_time = fig_time.add_subplot("111")
    ax_time.set_xlabel("Solar position (0 = periapsis) $[\pi RAD]$")
    ax_time.set_ylabel("Time $[s]$")
    ax_time.grid()

    theta = np.linspace(0, 2*pi, 1000)
    for e in ecc:
        orbit = CelestialOrbit(celestial_parent_body=Kerbin, inclination=0, longitude_ascending_node=0, argument_periapsis=0, eccentricity = e, semimajoraxis = 1e7)
        time = np.array([ orbit.TimeInShadow(t) for t in theta])
        ax_time.plot(((theta)  )/pi,time, label = "Ecc: " + str(e))

    ax_time.legend()
    plt.show()



    print("done")