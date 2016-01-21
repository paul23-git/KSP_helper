__author__ = 'Paul'
#from ksp_part_resource import *
import ksp_resources as KSP_RES
import ksp_part_resource as KSP_PART
import ksp_part_squad

class KSPModuleRTAntenna(KSP_RES.KSPModulePower):
    def __init__(self, power = 0, energy_cost = 0, omni_range=0, dish_range=0, **kwargs):
        super().__init__(name="ModuleRTAntenna", power_idle = energy_cost + power, **kwargs)
        self.__omni_range = omni_range
        self.__dish_range = dish_range

    def getAntennaType(self):
        return self.__dish_range == 0
    def isOmniAntenna(self):
        return self.__dish_range == 0
    def isDishAntenna(self):
        return self.__omni_range == 0

    @property
    def total_range(self):
        return max(self.__omni_range, self.__dish_range)

    @total_range.setter
    def total_range(self, value):
        if (self.isOmniAntenna()):
            self.__omni_range = value
        else:
            self.__dish_range = value

    @property
    def omni_range(self):
        return self.__omni_range

    @omni_range.setter
    def omni_range(self, value):
        self.__omni_range = value
        self.__dish_range = 0

    @property
    def dish_Range(self):
        return self.__dish_range

    @dish_Range.setter
    def dish_Range(self, value):
        self.__dish_range = value

    @property
    def energy_cost(self):
        return self.power_idle

    @energy_cost.setter
    def energy_cost(self, value):
        self.power_idle = value

def _load():
    _loadAntennas()
    _updateAntennas()

    print("RT update done")



def _loadAntennas():
    KSP_PART.Antenna(name="RTShortDish2",
            title="Reflectron KR-7",
            mass=0.5,
            price=800)
    KSP_PART.Antenna(name="RTShortAntenna1",
            title="Reflectron DP-10",
            mass=0.005,
            price=60)
    KSP_PART.Antenna(name="RTLongDish2",
            title="Reflectron KR-14",
            mass=1.0,
            price=2000)
    KSP_PART.Antenna(name="RTLongAntenna3",
            title="CommTech EXP-VR-2T",
            mass=0.02,
            price=400)
    KSP_PART.Antenna(name="RTLongAntenna2",
            title="Communotron 32",
            mass=0.01,
            price=600)
    KSP_PART.Antenna(name="RTGigaDish2",
            title="CommTech-1",
            mass=1,
            price=9500)
    KSP_PART.Antenna(name="RTGigaDish1",
            title="Reflectron GX-128",
            mass=0.5,
            price=11000)

def _updateAntennas():
    _updateRTAntennas()
    _updateSquadAntennas()

def _updateRTAntennas():
    try:
        p = KSP_PART.Antenna["RTShortAntenna1"]
        p.modify(tech_required="flightControl",
            antenna_module=KSPModuleRTAntenna(omni_range=500000.0, power=0.01))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTLongAntenna2"]
        p.modify(tech_required="largeElectrics",
            antenna_module=KSPModuleRTAntenna(omni_range=5000000.0, power=0.6))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTLongAntenna3"]
        p.modify(tech_required="specializedElectrics",
            antenna_module=KSPModuleRTAntenna(omni_range=3000000.0, power=0.18))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTShortDish1"]
        p.modify(antenna_module=KSPModuleRTAntenna(dish_range=90000000.0, power=0.82))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTShortDish2"]
        p.modify(tech_required="electrics",
            antenna_module=KSPModuleRTAntenna(omni_range=90000000.0, power=0.82))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTLongDish1"]
        p.modify(antenna_module=KSPModuleRTAntenna(omni_range=60000000000.0, power=0.93))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTLongDish2"]
        p.modify(tech_required="largeElectrics",
                 antenna_module=KSPModuleRTAntenna(omni_range=60000000000.0, power=0.93))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTGigaDish1"]
        p.modify(tech_required="flightControl",
            antenna_module=KSPModuleRTAntenna(omni_range=250000, power=0.13))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["RTGigaDish2"]
        p.modify(tech_required="flightControl",
            antenna_module=KSPModuleRTAntenna(omni_range=250000, power=0.13))
    except KeyError:
        pass

def _updateSquadAntennas():
    try:
        p = KSP_PART.Antenna["longAntenna"]
        p.modify(antenna_module=KSPModuleRTAntenna(omni_range=2500000, power=0.13))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["mediumDishAntenna"]
        p.modify(antenna_module=KSPModuleRTAntenna(dish_range=50000000, power=0.82))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["comDish"]
        p.modify(antenna_module=KSPModuleRTAntenna(dish_range=40000000000, power=0.93))
    except KeyError:
        pass
    try:
        p = KSP_PART.Antenna["HighGainAntenna"]
        p.modify(antenna_module=KSPModuleRTAntenna(dish_range=25000000000, power=1.04))
    except KeyError:
        pass

_load()