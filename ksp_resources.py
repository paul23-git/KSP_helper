import copy
import math as m


class KSPDict(dict):
    def __init__(self, seq = None, **kwargs):
        if seq is not None:
            if isinstance(seq, dict):
                super().__init__(seq)
                return
            else:
                try:
                    super().__init__(((p.name, p) for p in seq ))
                    return
                except TypeError:
                    pass
        super().__init__(**kwargs)

    def add(self, res):
        super().__setitem__(res.name, res)


KSPModuleKeywords = {
    "name",
    "power_idle",
    "power_active",
    "power",
    "nominalChargeRate",
    "chargeRate",
    "yawtorque",
    "pitchtorque",
    "rolltorque",
    "nominalChargeDistance",
    "sunTracking",
    "averageMultiplier"
}


class KSPModule(object):
    """
    KSP module base class
    """

    def __init__(self, name, **kwds):
        self.name = name
        super().__init__(**kwds)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

class KSPModulePower(KSPModule):
    """
    Module to for control operations
    """

    def __init__(self,  power_idle = 0, power_active = 0, **kwds):
        super().__init__( **kwds)
        self.__power_idle = power_idle
        self.__power_active = power_active
    @property
    def power_idle(self):
        return self.__power_idle

    @power_idle.setter
    def power_idle(self, t):
        self.__power_idle = t

    @property
    def power_active(self):
        return self.__power_active

    @power_active.setter
    def power_active(self, value):
        self.__power_active = value

    @property
    def power(self):
        return self.power_idle + self.power_active

    @power.setter
    def power(self, value):
        self.power_idle = value
        self.power_active = 0

class KSPModuleElectricityProducer(KSPModule):
    def __init__(self, nominalChargeRate, **kwargs):
        super().__init__(**kwargs)
        self.nominalChargeRate = nominalChargeRate

    @property
    def chargeRate(self):
        return self.nominalChargeRate

    @chargeRate.setter
    def chargeRate(self, value):
        self.nominalChargeRate = value

    def getChargeRate(self, distance=0):
        return self.chargeRate


class KSPModuleCommand(KSPModulePower):
    """
    Module to for control operations
    """

    def __init__(self, power, **kwds):
        super().__init__(name="ModuleCommand", power_idle=power, power_active=0, **kwds)


class KSPModuleReactionWheel(KSPModulePower):
    """
    Module for control operations
    """

    def __init__(self, power, yaw, pitch, roll, **kwds):
        super().__init__(power_idle=0, power_active=power, name = "ModuleReactionWheel", **kwds)
        self.yawtorque = yaw
        self.pitchtorque = pitch
        self.rolltorque = roll



class KSPModuleDeployableSolarPanel(KSPModuleElectricityProducer):
    """
    Module for solar panel action
    """
    def __init__(self, chargeRate, nominalChargeDistance=13599840256., sunTracking=True, **kwds):
        super().__init__(nominalChargeRate=chargeRate, name="ModuleDeployableSolarPanel", **kwds)
        self.sunTracking = sunTracking
        self.nominalChargeDistance = nominalChargeDistance
        self.averageMultiplier = 2/m.pi

    @property
    def chargeRate(self):
        v = self.nominalChargeRate*self.averageMultiplier
        if not self.sunTracking:
            v /= 2
        return v

    @chargeRate.setter
    def chargeRate(self, value):
        v = value/self.averageMultiplier
        if not self.sunTracking:
            v *= 2
        self.nominalChargeRate = v

    def getChargeRate(self, distance):
        return self.chargeRate * (self.nominalChargeDistance/distance)**2

class KSPModuleGenerator(KSPModuleElectricityProducer):
    def __init__(self, chargeRate, **kwargs):
        super().__init__(nominalChargeRate=chargeRate, name="ModuleGenerator", **kwargs)

class KSPResource(object):
    """
    Generig resource class
    """
    def __init__(self, name, density=0, unit_cost=0, specific_heat=0):
        """

        :type name: str
        :param name: id
        :param density: kg/unit
        :type density: float
        :param unit_cost: cost/unit
        :type unit_cost: float
        :param specific_heat: specific heat
        :type specific_heat: float
        """
        self.__name = name
        self.density = density
        self.unit_cost = unit_cost
        self.specific_heat = specific_heat
    @property
    def name(self):
        return self.__name

    def __hash__(self):
        return hash(self.__name)

    def __eq__(self, other):
        return self.__name == other.name

ALL_RESOURCES = KSPDict({
    KSPResource("ElectricCharge"),
    KSPResource("LiquidFuel", 0.005, 0.8, 2100),
    KSPResource("Oxidizer", 0.005, 0.18, 1551),
    KSPResource("SolidFuel", 0.0075, 0.6, 920),
    KSPResource("MonoPropellant", 0.004, 1.2, 3000),
    KSPResource("XenonGass", 0.0001, 4, 120),
    KSPResource("Ore", 0.01, 0.02, 1000),
    KSPResource("Ablator", 0.001, 0.5, 400)
})

__all__ = [
    "KSPModule",
    "KSPModulePower",
    "KSPModuleCommand",
    "KSPModuleElectricityProducer",
    "KSPModuleDeployableSolarPanel",
    "KSPModuleGenerator",
    "KSPModuleReactionWheel",
    "ALL_RESOURCES"
]