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
import copy
import itertools
import math

class satellite(object):
    def __init__(self, orbit, payload=None, structural_mass=0, structural_price=0):
        self.__all_parts = []
        self.orbit = orbit
        self.payload = payload
        self.structural_mass = structural_mass
        self.structural_price = structural_price

    #def __getattr__(self, item):
    #    found = False
    #    v = 0
    #    for i in self.__all_parts:
    #        t = getattr(i, item, None)
    #        if t is not None:
    #            v += t
    #            found = True
    #    if found:
    #        return v
    #    else:
    #        raise AttributeError(item)

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

def _add_batteries_to_sat(satellite, base_battery, energy_storage_reuirements):
    n = max(0, math.floor(energy_storage_reuirements/base_battery.electric_charge))
    for i in range(n):
        satellite.addPart(base_battery)
    return n * base_battery.electric_charge

def _generate_sat_with_batteries(base_satellite, battery, energy_storage_requirements, possible_battery_list, first_order_check):
    print(base_satellite)
    newsat = copy.copy(base_satellite)
    energy_storage_requirements -= _add_batteries_to_sat(newsat, battery, energy_storage_requirements)
    if energy_storage_requirements > 0:
        sats = []
        possible_battery_list = [b for b in possible_battery_list if b.electric_charge <= battery.electric_charge]
        bad_storage = set(c[0] for c in itertools.permutations(possible_battery_list, 2) if not first_order_check(*c))
        useful_storage = (s for s in possible_battery_list if s not in bad_storage)
        for b in useful_storage:
            sats.extend(_generate_sat_with_batteries(newsat, b, energy_storage_requirements, possible_battery_list, first_order_check))
        return sats
    else:
        return [newsat]




def generate_all_battery_sattelites(base_satellite, energy_storage_requirements, tech_level, first_order_check = None):
    storage_list = [battery for battery in KSP_PART.BatteryPack if battery.tech_required in tech_level]
    storage_list.sort(key=lambda p: p.title,reverse=False)
    print([s.title for  s in storage_list])
    if first_order_check is None:
        useful_storage = storage_list
    else:
        #check for bad elements
        bad_storage = set(c[0] for c in itertools.permutations(storage_list, 2) if not first_order_check(*c) )
        useful_storage = [s for s in storage_list if s not in bad_storage]
    for battery in useful_storage:
        print(battery.title)

    print(useful_storage)
    all_copies = itertools.chain([_generate_sat_with_batteries(base_satellite,
                                                             battery,
                                                             energy_storage_requirements,
                                                             storage_list,
                                                             first_order_check)
                  for battery in useful_storage])



