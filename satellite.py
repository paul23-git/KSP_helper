__author__ = 'Paul'

import ksp_resources
import numpy as np
import math as m
import scipy as sp
import matplotlib.pyplot as plt
import ksp_part_resource as KSP_PART
import ksp_part_squad
import ksp_part_remotetech
import celestial_body
import functools
import operator


class satellite(object):
    def __init__(self, orbit, payload=None, structural_mass=0, structural_price=0):
        self.__all_parts = []
        self.orbit = orbit
        self.payload = payload
        self.structural_mass = structural_mass
        self.structural_price = structural_price

    def __getattr__(self, item):
        return sum(getattr(i, item, 0) for i in self.__all_parts)

    def addPart(self, part):
        self.__all_parts.append(part)

    def getTotalPowerConsumption(self):
        p = sum(i.power for i in self.__all_parts)
        return p

    def getIdlePowerConsumption(self):
        return sum(i.power_idle for i in self.__all_parts)

    def getTotalWeight(self):
        return sum(i.mass for i in self.__all_parts)

    def getTotalEnergyStorage(self):
        return sum(i.energy for i in self.__all_parts)

    def getPrice(self):
        return sum(i.price for i in self.__all_parts)

    def getChargeRate(self, distance):
        return sum(i.getChargeRate(distance) for i in self.__all_parts)

    def getAverageChargeRate(self):
        d = self.getMaxAverageDistance()
        return self.getChargeRate(d)

    def getMinChargeRate(self):
        d = self.getMaxDistance()
        return self.getChargeRate(d)

    def getMaxAverageDistance(self):
        o = self.orbit
        while not isinstance(o.parent, celestial_body.CelestialStar):
            o = o.parent.orbit
        return o.apoapsis_distance

    def getMaxDistance(self):
        o = self.orbit
        d = o.apoapsis_distance
        while not isinstance(o.parent, celestial_body.CelestialStar):
            o = o.parent.orbit
            d += o.apoapsis_distance
        return d

    @property
    def price(self):
        price = self.structural_price + sum(part.price for part in self.__all_parts if hasattr(part, "price"))
        try:
            price += self.payload.price
        except AttributeError:
            pass
        return price

    @property
    def mass(self):
        mass = self.empty_mass
        try:
            mass += self.payload.mass
        except AttributeError:
            pass
        return mass

    @property
    def empty_mass(self):
        return self.structural_mass + sum(part.mass for part in self.__all_parts if hasattr(part, "mass"))


def calculate_minimal_electric_storage(satellite):
     return satellite.getIdlePowerConsumption() * satellite.orbit.TotalMaxTimeInShadow()

def calculate_minimal_charge_rate(satellite):
    P = satellite.getIdlePowerConsumption()
    for t in satellite.orbit.GenerateShadowLightTimeList():
        P *= (t[0] + t[1]) / t[0]
    return P