__author__ = 'Paul'

import copy
from ksp_resources import ALL_RESOURCES, KSPModule, KSPModuleKeywords

class IterPartRegistry(type):
    def __init__(cls, name, bases, attrs):
        super().__init__(name, bases, attrs)
        cls._registry = {}

    def __contains__(cls, item):
        return (cls._registry.__contains__(item) or
                any(c.__contains__(item) for c in cls.__subclasses__()))

    def __iter__(cls):
        # return (c for c in cls._registry.values() if isinstance(c, cls))
        yield from cls._registry.values()
        for subcls in cls.__subclasses__():
            yield from subcls

    def __call__(cls, name, *args, **kwargs):
        if name not in KSPPart:
            ret = super().__call__(name=name, *args, **kwargs)
            cls._registry[name] = ret
        else:
            ret = KSPPart[name]
            if isinstance(ret, cls):
                ret.modify(**kwargs)
            else:
                raise TypeError(name + " is " +str(type(ret)) + ", expected " + str(cls))
        return ret

    def __getitem__(cls, item):
        try:
            return cls._registry[item]
        except KeyError:
            for subcls in cls.__subclasses__():
                try:
                    return subcls.__getitem__(item)
                except KeyError:
                    pass
            raise

    def __len__(cls):
        return len(cls._registry) + sum(c.__len__() for c in cls.__subclasses__())

    def __delitem__(cls, key):
        try:
            del cls._registry[key]
            return
        except KeyError:
            for subcls in cls.__subclasses__():
                try:
                    subcls.__delitem__(key)
                    return
                except KeyError:
                    pass
            raise


class KSPPart(object, metaclass=IterPartRegistry):
    """
    Main part class
    """

    _resourceNameBindings = {}

    def __getattr__(self, item):
        if item in self._resourceNameBindings:
            return self.resources[self._resourceNameBindings[item]]
        if item in KSPModuleKeywords:
            return sum(getattr(m, item, 0) for m in self.modules)
        raise AttributeError(item)

    def __init__(self, name, mass, techrequired="start", title=None, price=0, extra_modules=None, resources=None, enabled=True):
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
        self.__name = name
        # self._registry.append(self)
        self.mass = mass
        self.price = price
        if title is None:
            self.title = name
        else:
            self.title = title
        if extra_modules is not None:
            self.modules = copy.deepcopy(extra_modules)
        else:
            self.modules = set()
        if resources is not None:
            self.resources = copy.copy(resources)
        else:
            self.resources = {}
        self.enabled = enabled
        self.tech_required = techrequired

    def modify(self, **kwargs):
        for k, v in kwargs.items():
            if k in self._resourceNameBindings:
                self.resources[self._resourceNameBindings[k]] = v
            elif isinstance(v, KSPModule):
               self.modules.discard(v)
               self.modules.add(v)
            else:
                setattr(self, k, v)

    def hasResource(self, name):
        if name in self._resourceNameBindings:
            name = self._resourceNameBindings[name]
        return name in self.resources

    def getResource(self, name):
        if name in self._resourceNameBindings:
            name = self._resourceNameBindings[name]
        return self.resources[name]

    def hasModule(self, name):
        return any(m.name == name for m in self.modules)

    def getModule(self, name):
        return next(m for m in self.modules if m.name == name)

    @property
    def name(self):
        return self.__name

    @property
    def energy(self):
        try:
            return self.resources["ElectricCharge"]
        except KeyError:
            return 0

    @property
    def total_mass(self):
        return self.mass + sum(ALL_RESOURCES[k].density * v for k, v in self.resources.items())

    @property
    def total_price(self):
        return self.price + sum(ALL_RESOURCES[k].price * v for k, v in self.resources.items())

    def __hash__(self):
        return hash((self.name, type(self)))

    def __eq__(self, other):
        return type(self) == type(other) and self.name == other.name

    def __str__(self, *args, **kwargs):
        return self.name


class StructuralPart(KSPPart):
    pass

class PayloadPart(KSPPart):
    pass

class ControlPart(KSPPart):
    """
    Control probe class
    """
    _resourceNameBindings = copy.copy(KSPPart._resourceNameBindings)
    _resourceNameBindings.update(
        {"ElectricCharge":"ElectricCharge"}
    )
    def __init__(self, name, mass, command_module, controlwheel_module=None, ElectricCharge=0, crew_capacity=0, **kwds):
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
        :param ElectricCharge: Amount of energy
        :type ElectricCharge: float
        """
        super().__init__(name=name, mass=mass, **kwds)
        self.modules.add(copy.deepcopy(command_module))
        self.resources["ElectricCharge"] = ElectricCharge
        if controlwheel_module is not None:
            self.modules.add(copy.deepcopy(controlwheel_module))
        self.crew_capacity = crew_capacity


class ControlProbe(ControlPart):
    pass


class ControlCommand(ControlPart):
    _resourceNameBindings = copy.copy(ControlPart._resourceNameBindings)
    _resourceNameBindings.update(
        {"mono_prop_amount":"Mono"}
    )
    def __init__(self, name, mass, command_module, controlwheel_module, ElectricCharge=0, mono_prop_amount=0,
                 crew_capacity=0, **kwds):
        super().__init__(name=name, mass=mass, command_module=command_module, controlwhell_module=controlwheel_module,
                         ElectricCharge=ElectricCharge, crew_capacity=crew_capacity, **kwds)
        self.resources["Mono"] = mono_prop_amount


class BatteryPack(KSPPart):
    """
    Control probe class
    """
    _resourceNameBindings = copy.copy(KSPPart._resourceNameBindings)
    _resourceNameBindings.update(
        {"electric_charge":"ElectricCharge"}
    )
    def __init__(self, name, mass, electric_charge, **kwds):
        super().__init__(name=name, mass=mass, **kwds)
        self.resources["ElectricCharge"] = electric_charge


class EnergyGenerationPart(KSPPart):
    def __init__(self, name, mass, ModuleEnergy, **kwargs):
        super().__init__(name=name, mass=mass, **kwargs)
        self.modules.add(copy.deepcopy(ModuleEnergy))

    def getChargeRate(self, distance=0):
        return sum(m.getChargeRate(distance) for m in self.modules if hasattr(m,"getChargeRate"))


class Antenna(KSPPart):
    pass


def ModifyOrAdd(cls, name, pred=None, **kwargs):
    if pred is None:
        pred = lambda x: x.name == name

    all_old = (x for x in cls if pred(x))
    done = False
    for old in all_old:
        done = True
        old.modify(**kwargs)
    if done:
        return None
    else:
        return cls(name=name, **kwargs)


def Modify(pred=None, name=None, **kwargs):
    if pred is None:
        if name is not None:
            pred = lambda x: x.name == name
        else:
            pred = lambda x: True
    try:
        all_old = (old for old in KSPPart if pred(old))
        for old in all_old:
            old.modify(**kwargs)
    except StopIteration:
        return None
