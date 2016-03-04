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


class satellite(celestial_body.CelestialBody_Orbitting):
    def __init__(self, orbit, payload=None, structural_mass=0, structural_price=0, *args, **kwargs):
        super(mass=0, mu=0, radius=0, sidereal_rotation_period=np.inf, SOI_radius=0, orbit=orbit, *args, **kwargs)
        self._all_parts = []
        self.payload = payload
        self.structural_mass = structural_mass
        self.structural_price = structural_price

    def __getattr__(self, item):
        if item.startswith('_'):
            raise AttributeError(item)
        found = False
        v = 0
        for i in self._all_parts:
            t = getattr(i, item, None)
            if t is not None:
                v += t
                found = True
        if found:
            return v
        else:
            raise AttributeError(item)

    def __copy__(self):
        s = type(self)(self.orbit, self.payload, self.structural_mass, self.structural_price)
        s._all_parts = self.getCopyOfPartlist()
        return s

    def copy(self):
        return self.__copy__()

    def addPart(self, part):
        self._all_parts.append(part)

    def getNumberParts(self, type=None):
        if type is None:
            return len(self._all_parts)
        else:
            n = 0
            for _ in self._all_parts:
                n += isinstance(_,type)
            return n

    def getCopyOfPartlist(self, pred=None, *args, **kwargs):
        if pred is None:
            return copy.copy(self._all_parts)
        else:
            return [p for p in self._all_parts if pred(p, *args, **kwargs)]

    def getTotalPowerConsumption(self):
        p = sum(i.power for i in self._all_parts)
        return p

    def getIdlePowerConsumption(self):
        return sum(i.power_idle for i in self._all_parts)

    def getTotalWeight(self):
        return sum(i.mass for i in self._all_parts)

    def getTotalEnergyStorage(self):
        return sum(i.energy for i in self._all_parts)

    def getPrice(self):
        return sum(i.price for i in self._all_parts)

    def getChargeRate(self, distance):
        return sum(i.getChargeRate(distance) for i in self._all_parts)

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
        price = self.structural_price + sum(part.price for part in self._all_parts if hasattr(part, "price"))
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
        return self.structural_mass + sum(part.mass for part in self._all_parts if hasattr(part, "mass"))


def calculate_minimal_electric_storage(satellite):
    return satellite.getIdlePowerConsumption() * satellite.orbit.total_max_time_in_shadow()


def calculate_minimal_charge_rate(satellite):
    P = satellite.getIdlePowerConsumption()
    for t in satellite.orbit.generate_shadow_light_time_list():
        P *= (t[0] + t[1]) / t[0]
    return P


def _add_batteries(satellite, base_battery, energy_storage_requirements):
    n = max(1, math.floor(energy_storage_requirements / base_battery.electric_charge))
    for i in range(n):
        satellite.addPart(base_battery)
    return n


def _generate_batteries(satellite, energy_storage_requirements, possible_battery_list, first_order_check,
                        battery=None):
    if battery is not None:
        satellite = copy.copy(satellite)
        energy_storage_requirements -= _add_batteries(satellite, battery,
                                                      energy_storage_requirements) * battery.electric_charge

    if energy_storage_requirements > 0:
        sats = []
        if battery is not None:
            possible_battery_list = [b for b in possible_battery_list if b.electric_charge <= battery.electric_charge]
        if first_order_check is None:
            useful_storage = possible_battery_list
        else:
            useful_storage = _test_func2(possible_battery_list, first_order_check, energy_storage_requirements)
        for b in useful_storage:
            newsat = _generate_batteries(satellite, energy_storage_requirements, possible_battery_list,
                                         first_order_check, battery=b)
            sats.extend(newsat)
        return sats
    else:
        return [satellite]


def _test_func(all_copies, func, *args, **kwargs):
    bad_copies = set(c[0] for c in itertools.permutations(all_copies, 2) if not func(c[0], c[1], *args, **kwargs))
    all_copies = (c for c in all_copies if c not in bad_copies)
    return all_copies
def _test_func2(all_copies, func, *args, **kwargs):
    for i in all_copies:
        if all(i == j or func(i,j, *args, **kwargs) for j in all_copies):
            yield i

def generate_all_battery(base_satellite, tech_level, energy_storage_requirements = None, first_order_check=None):
    if energy_storage_requirements is None:
        energy_storage_requirements = calculate_minimal_electric_storage(base_satellite)
    storage_list = [battery for battery in KSP_PART.BatteryPack if battery.tech_required in tech_level]
    storage_list.sort(key=lambda p: p.title, reverse=True)
    all_copies = _generate_batteries(base_satellite,
                                     energy_storage_requirements,
                                     storage_list,
                                     first_order_check)

    def check_post(sat_l, sat_r):
        left_v = (sat_l.mass, sat_l.price, sat_l.getNumberParts())
        right_v = (sat_r.mass, sat_r.price, sat_r.getNumberParts())
        return any(ll < rr for ll, rr in zip(left_v, right_v))

    all_copies = list(_test_func2(all_copies, check_post))
    return all_copies


def _add_solar_cells(satellite, solar_cell, power_requirements, distance):
    n = max(1, math.floor(power_requirements / solar_cell.getChargeRate(distance)))
    for i in range(n):
        satellite.addPart(solar_cell)
    return n


def _generate_solar_cells(satellite, power_requirements, possible_cell_list,
                          first_order_check, distance, solar_cell=None, *args, **kwargs):
    if solar_cell is not None:
        satellite = copy.copy(satellite)
        power_requirements -= _add_solar_cells(satellite, solar_cell, power_requirements,
                                               distance) * solar_cell.getChargeRate(distance)

    if power_requirements:
        sats = []
        if solar_cell is not None:
            oldRate = solar_cell.getChargeRate(distance)
            possible_cell_list = [cell for cell in possible_cell_list if cell.getChargeRate(distance) >= oldRate]
        if first_order_check is None:
            useful_cells = possible_cell_list
        else:
            useful_cells = _test_func(possible_cell_list, first_order_check,
                                      power_requirements, distance, *args, **kwargs)
        for cell in useful_cells:
            newsat = _generate_solar_cells(satellite, power_requirements, possible_cell_list,
                                           first_order_check, distance, solar_cell=cell)
            sats.extend(newsat)
        return sats
    else:
        return [satellite]


def generate_all_solar_cells(base_satellite, tech_level,  power_requirements = None,
                             distance=None, first_order_check=None, *args, **kwargs):
    if power_requirements is None:
        power_requirements = calculate_minimal_charge_rate(base_satellite)
    if distance is None:
        distance = base_satellite.getAvgDistance()
    solar_cell_list = [solar_cell for solar_cell in KSP_PART.EnergyGenerationPart if
                       solar_cell.tech_required in tech_level]
    all_copies = _generate_solar_cells(base_satellite,
                                       power_requirements,
                                       solar_cell_list,
                                       first_order_check,
                                       distance,
                                       *args, **kwargs)
    return all_copies
