import copy


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

class IterPartRegistry(type):
    def __iter__(cls):
        return iter(cls._registry)
class KSPPart(object):
    """
    Main part class
    """

    __metaclass__ = IterPartRegistry
    _registry = []

    def __init__(self, name, mass, cost=0, extra_modules=None, resources=None):
        """
        Creates basic part
        :param name: id
        :type name: str
        :param mass: mass
        :type mass: float
        :param extra_modules: set containing all extra modules.
        :param resources: dictionary containing resource, amount keys
        :type resources: {KSPResource: float}
        """
        self.name = name
        self._registry.append(self)
        self.mass = mass
        if extra_modules is not None:
            self.modules = copy.deepcopy(extra_modules)
        else:
            self.modules = set()
        if resources is not None:
            self.resources = copy.copy(resources)
        else:
            self.resources = {}
    @property
    def energy(self):
        if "ElectricalCharge" in self.resources:
            v = self.resources["ElectricalCharge"]
            return v


    def __hash__(self):
        return hash((self.name, type(self)))

    def __eq__(self, other):
        return type(self) == type(other) and self.name == other.name

class KSPModule(object):
    """
    KSP module base class
    """

    def __init__(self, name):
        self.name = name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

class ModuleCommand(KSPModule):
    """
    Module to for control operations
    """

    def __init__(self, power):
        super().__init__("ModuleCommand")
        self.power = power

class ModuleReactionWheel(KSPModule):
    """
    Module for control operations
    """

    def __init__(self, power, yaw, pitch, roll):
        super().__init__("ModuleReactionWheel")
        self.power = power
        self.yawtorque = yaw
        self.pitchtorque = pitch
        self.rolltorque = roll

class ModuleDeployableSolarPanel:
    """
    Module for solar panel action
    """
    def __init__(self, chargeRate, sunTracking=True):
        self.chargeRate = chargeRate
        self.sunTracking = sunTracking

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

class ControlProbe(KSPPart):
    """
    Control probe class
    """

    def __init__(self, name, mass, command_module, controlwheel_module, charge_amount):
        """
        Creats control probe
        :param name: ID
        :type name: str
        :param mass: mass
        :type mass: float
        :param command_module: Command module (copied)
        :type command_module: ModuleCommand
        :param controlwheel_module: sas module (copyied)
        :type controlwheel_module: ModuleReactionWheel
        :param charge_amount: Amount of energy
        :type charge_amount: float
        """
        super().__init__(name, mass)
        self.modules.add(copy.deepcopy(command_module))
        r = ALL_RESOURCES["ElectricCharge"]
        self.resources[r] =  charge_amount
        self.modules.add(copy.deepcopy(controlwheel_module))


class BatteryPack(KSPPart):
    """
    Control probe class
    """

    def __init__(self, name, mass, charge_amount):
        super().__init__(name, mass)
        self.resources[ALL_RESOURCES["ElectricCharge"]] =  charge_amount
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