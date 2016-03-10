__author__ = 'Paul'

import ksp_resources
import numpy as np
import math as m
import scipy as sp
import matplotlib.pyplot as plt
from ksp_part_resource import KSPPart, KSPModule, KSPModuleKeywords, BatteryPack, EnergyGenerationPart
import ksp_part_squad
import ksp_part_remotetech
import celestial_body as CB
import celestial_orbit as CO
import functools
import operator
import copy
import itertools
import math
import typing
from typing import Callable, List, Iterable, TypeVar, Generator, Container, Sequence, MutableSequence
T = TypeVar('T')
TPredFunc = Callable[...,bool]

class Satellite(CB.CelestialObjectAbstract_Orbitting):
    def __init__(self, orbit, payload=None, structural_mass=0, structural_price=0, *args, **kwargs):
        super().__init__(orbit=orbit, *args, **kwargs)
        self._all_parts = [] #type:List[KSPPart]
        self.payload = payload #type:float
        self.structural_mass = structural_mass #type:float
        self.structural_price = structural_price #type:float

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
        s = type(self)(orbit=self.orbit, payload=self.payload, structural_mass=self.structural_mass,
                       structural_price=self.structural_price, name=self.name)
        s._all_parts = self.getCopyOfPartlist()
        return s

    def copy(self):
        return self.__copy__()

    def addPart(self, part:KSPPart):
        self._all_parts.append(part)

    def getNumberParts(self, type:KSPPart=None) -> float:
        if type is None:
            return len(self._all_parts)
        else:
            n = 0
            for _ in self._all_parts:
                n += isinstance(_,type)
            return n

    def getCopyOfPartlist(self, pred:Callable[...,bool]=None, *args, **kwargs):
        if pred is None:
            return copy.copy(self._all_parts)
        else:
            return [p for p in self._all_parts if pred(p, *args, **kwargs)]

    def getTotalPowerConsumption(self) -> float:
        p = sum(i.power for i in self._all_parts)
        return p

    def getIdlePowerConsumption(self) -> float:
        return sum(i.power_idle for i in self._all_parts)

    def getTotalWeight(self) -> float:
        return sum(i.mass for i in self._all_parts)

    def getTotalEnergyStorage(self) -> float:
        return sum(i.energy for i in self._all_parts)

    def getPrice(self) -> float:
        return sum(i.price for i in self._all_parts)

    def getChargeRate(self, distance) -> float:
        return sum(i.getChargeRate(distance) for i in self._all_parts)


    def getMinChargeRate(self) -> float:
        d = self.max_distance_to_star()
        return self.getChargeRate(d)

    def average_max_distance_to_star(self):
        o = self.orbit
        while not isinstance(o.parent, CB.CelestialStar):
            o = o.parent.orbit
        return o.apoapsis_distance

    def max_distance_to_star(self, eps:float=0.00001) -> float:
        p = self.parent #type: CB.CelestialBody
        while not isinstance(p, CB.CelestialStar):
            p = p.parent
        d = self.orbit.get_total_max_distance(p, eps=eps)
        return d

    @property
    def price(self) -> float:
        price = self.structural_price + sum(part.price for part in self._all_parts if hasattr(part, "price"))
        try:
            price += self.payload.price
        except AttributeError:
            pass
        return price

    @property
    def mass(self) -> float:
        mass = self.empty_mass
        try:
            mass += self.payload.mass
        except AttributeError:
            pass
        return mass

    @property
    def empty_mass(self) -> float:
        return self.structural_mass + sum(part.mass for part in self._all_parts if hasattr(part, "mass"))


    def calculate_minimal_electric_storage(self) -> float:
        return self.getIdlePowerConsumption() * self.orbit.iterative_max_time_in_shadow()


def calculate_minimal_charge_rate(satellite:Satellite) -> float:
    P = satellite.getIdlePowerConsumption()
    for t in satellite.orbit.generate_shadow_light_time_list():
        P *= (t[0] + t[1]) / t[0]
    return P


def _add_batteries(battery_list:MutableSequence[BatteryPack], base_battery:BatteryPack,
                   energy_storage_requirements:float) -> int:
    n = max(1, math.floor(energy_storage_requirements / base_battery.electric_charge))
    for i in range(n):
        battery_list.append(base_battery)
    return n


def _generate_batteries(already_added_batteries:MutableSequence[BatteryPack],
                        energy_storage_requirements:float, possible_battery_list:Iterable[BatteryPack],
                        first_order_check:TPredFunc=None, battery:BatteryPack=None) -> List[BatteryPack]:
    if battery is not None:
        already_added_batteries = copy.copy(already_added_batteries)
        energy_storage_requirements -= _add_batteries(already_added_batteries, battery,
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
            newlist = _generate_batteries(already_added_batteries, energy_storage_requirements,
                                         possible_battery_list, first_order_check, battery=b)
            sats.extend(newlist)
        return sats
    else:
        return [already_added_batteries]

def _test_func2(all_copies:Iterable[T], func:TPredFunc, *args, **kwargs) -> Generator[T, None, None]:
    for i in all_copies:
        if all(i == j or func(i,j, *args, **kwargs) for j in all_copies):
            yield i


def generate_all_battery(tech_level:Container[str],
                         energy_storage_requirements:float, first_order_check:TPredFunc=None) -> List[Satellite]:
    storage_list = [battery for battery in BatteryPack if battery.tech_required in tech_level]
    storage_list.sort(key=lambda p: p.title, reverse=True)
    all_copies = _generate_batteries([],
                                     energy_storage_requirements,
                                     storage_list,
                                     first_order_check)

    def check_post(add_l:Sequence[BatteryPack], add_r:Sequence[BatteryPack], eps:float=0.0000001) -> bool:
        t = tuple((i.mass, i.price) for i in add_l)
        left_v = [sum(x) for x in zip(*t)] + [len(add_l)]
        t = tuple((i.mass, i.price) for i in add_r)
        right_v = [sum(x) for x in zip(*t)] + [len(add_r)]
        return any((rr - ll) > ll*eps  for ll, rr in zip(left_v, right_v))
    #KSP-SPECIFIC: remove too many parts items
    min_parts = min(len(i) for i in all_copies)
    all_copies = [copy for copy in all_copies if len(copy) <= 10*min_parts]

    all_copies = list(_test_func2(all_copies, check_post))
    return all_copies


def _add_solar_cells(cell_list:MutableSequence[EnergyGenerationPart], solar_cell:EnergyGenerationPart,
                     power_requirements:float, distance:float) -> int:
    n = max(1, math.floor(power_requirements / solar_cell.getChargeRate(distance)))
    for i in range(n):
        cell_list.append(solar_cell)
    return n


def _generate_solar_cells(already_added_cells:MutableSequence[EnergyGenerationPart],
                          power_requirements:float, possible_cell_list:Iterable[EnergyGenerationPart],
                          first_order_check:TPredFunc, distance:float, solar_cell:EnergyGenerationPart=None,
                          *args, **kwargs):
    if solar_cell is not None:
        already_added_cells = copy.copy(already_added_cells)
        power_requirements -= _add_solar_cells(already_added_cells, solar_cell, power_requirements,
                                               distance) * solar_cell.getChargeRate(distance)

    if power_requirements > 0:
        sats = []
        if solar_cell is not None:
            oldRate = solar_cell.getChargeRate(distance)
            possible_cell_list = [cell for cell in possible_cell_list if cell.getChargeRate(distance) >= oldRate]
        if first_order_check is None:
            useful_cells = possible_cell_list
        else:
            useful_cells = _test_func2(possible_cell_list, first_order_check,
                                      distance, power_requirements, *args, **kwargs)
        for cell in useful_cells:
            newlist = _generate_solar_cells(already_added_cells, power_requirements, possible_cell_list,
                                           first_order_check, distance, solar_cell=cell)
            sats.extend(newlist)
        return sats
    else:
        return [already_added_cells]


def generate_all_solar_cells(tech_level:Container[str], power_requirements:float,
                             distance:float, first_order_check:TPredFunc=None, *args, **kwargs):
    solar_cell_list = [solar_cell for solar_cell in EnergyGenerationPart if
                       solar_cell.tech_required in tech_level]
    for j in solar_cell_list:
        print(j.title)
    all_copies = _generate_solar_cells([],
                                       power_requirements,
                                       solar_cell_list,
                                       first_order_check,
                                       distance,
                                       *args, **kwargs)
    def check_post(add_l:Sequence[EnergyGenerationPart],
                   add_r:Sequence[EnergyGenerationPart], eps:float=0.0000001) -> bool:
        t = tuple((i.mass, i.price) for i in add_l)
        left_v = [sum(x) for x in zip(*t)] + [len(add_l)]
        t = tuple((i.mass, i.price) for i in add_r)
        right_v = [sum(x) for x in zip(*t)] + [len(add_r)]
        print(left_v, right_v, [(rr - ll) > ll*eps for ll, rr in zip(left_v, right_v)])
        return any((rr - ll) > ll*eps  for ll, rr in zip(left_v, right_v))

    #KSP-SPECIFIC: remove too many parts items
    min_parts = min(len(i) for i in all_copies)
    all_copies = [copy for copy in all_copies if len(copy) <= 10*min_parts]

    all_copies = list(_test_func2(all_copies, check_post))
    return all_copies
