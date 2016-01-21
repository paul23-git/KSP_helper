__author__ = 'Paul'
from ksp_part_resource import *
from ksp_resources import KSPModuleCommand, \
    KSPModuleReactionWheel, \
    KSPModuleDeployableSolarPanel, \
    KSPModuleGenerator


def _load():
    _loadEnergy()
    _loadCommand()
    _loadAntennas()

def _loadEnergy():
    EnergyGenerationPart(name="solarPanels5",
                         title="OX-STAT Photovoltaic Panels",
                         mass=0.005,
                         price=75,
                         ModuleEnergy=KSPModuleDeployableSolarPanel(chargeRate=0.35,
                                                                 sunTracking=False),
                         techrequired="electrics")

    EnergyGenerationPart(name="solarPanels4",
                         title="OX-4L 1x6 Photovoltaic Panels",
                         mass=0.0175,
                         price=380,
                         ModuleEnergy=KSPModuleDeployableSolarPanel(chargeRate=1.64),
                         techrequired="advElectrics")

    EnergyGenerationPart(name="largeSolarPanel",
                         title="Gigantor XL Solar Array",
                         mass=0.3,
                         price=3000,
                         ModuleEnergy=KSPModuleDeployableSolarPanel(chargeRate=24.4),
                         techrequired="largeElectrics")

    EnergyGenerationPart(name="rtg",
                         title="PB-NUK Radioisotope Thermoelectric Generator",
                         mass=0.08,
                         price=23300,
                         ModuleEnergy=KSPModuleGenerator(chargeRate=0.75),
                         techrequired="experimentalElectrics")

    BatteryPack(name="batteryPack",
                title="Z-100 Rechargeable Battery Pack",
                mass=0.005,
                price=80,
                ElectricCharge=100,
                techrequired="basicScience")

    BatteryPack(name="batteryBankMini",
                title="Z-200 Rechargeable Battery Pack",
                mass=0.01,
                price=360,
                ElectricCharge=200,
                techrequired="electrics")

    BatteryPack(name="ksp_r_largeBatteryPack",
                title="Z-400 Rechargeable Battery Pack",
                mass=0.02,
                price=550,
                ElectricCharge=400,
                techrequired="advElectrics")

    BatteryPack(name="batteryBank",
                title="Z-1k Rechargeable Battery Pack",
                mass=0.05,
                price=880,
                ElectricCharge=1000,
                techrequired="largeElectrics")

    BatteryPack(name="batteryBankLarge",
                title="Z-4k Rechargeable Battery Pack",
                mass=0.2,
                price=4500,
                ElectricCharge=4000,
                techrequired="specializedElectrics")


def _loadCommand():
    ControlProbe(name="probeCoreSphere",
                 title="Stayputnik Mk. 1",
                 mass=0.05,
                 price=300,
                 command_module=KSPModuleCommand(power = 0.0277778),
                 ElectricCharge=10,
                 techrequired="basicScience")
    ControlProbe(name="probeStackSmall",
                 title="RC-001S Remote Guidance Unit",
                 mass=0.1,
                 price=2250,
                 command_module=KSPModuleCommand(power=0.05),
                 controlwheel_module=KSPModuleReactionWheel(power=0.03,
                                                         yaw=0.5,
                                                         pitch=0.5,
                                                         roll=0.5),
                 ElectricCharge=15,
                 techrequired="advUnmanned")
    ControlProbe(name="probeStackLarge",
                 title="RC-L01 Remote Guidance Unit",
                 mass=0.5,
                 price=3400,
                 command_module=KSPModuleCommand(power=0.08),
                 controlwheel_module=KSPModuleReactionWheel(power=0.15,
                                                         yaw=1.5,
                                                         pitch=1.5,
                                                         roll=1.5),
                 ElectricCharge=30,
                 techrequired="largeUnmanned")
    ControlProbe(name="roverBody",
                 title="Probodobodyne RoveMate",
                 mass=0.15,
                 price=800,
                 command_module=KSPModuleCommand(power=0.04),
                 ElectricCharge=120,
                 techrequired="largeUnmanned")
    ControlProbe(name="probeCoreOcto",
                 title="Probodobodyne OKTO",
                 mass=0.1,
                 price=450,
                 command_module=KSPModuleCommand(power=0.02),
                 controlwheel_module=KSPModuleReactionWheel(power=0.03,
                                                         yaw=0.3,
                                                         pitch=0.3,
                                                         roll=0.3),
                 ElectricCharge=10,
                 techrequired="electrics")
    ControlProbe(name="probeCoreOcto2",
                 title="Probodobodyne OKTO2",
                 mass=0.04,
                 price=1480,
                 command_module=KSPModuleCommand(power=0.03),
                 ElectricCharge=5,
                 techrequired="unmannedTech")
    ControlProbe(name="probeCoreHex",
                 title="Probodobodyne HECS",
                 mass=0.1,
                 price=650,
                 command_module=KSPModuleCommand(power=0.025),
                 controlwheel_module=KSPModuleReactionWheel(power=0.05,
                                                         yaw=0.5,
                                                         pitch=0.5,
                                                         roll=0.5),
                 ElectricCharge=10,
                 techrequired="precisionEngineering")
    ControlProbe(name="probeCoreCub",
                 title="Probodobodyne QBE",
                 mass=0.07,
                 price=360,
                 command_module=KSPModuleCommand(power=0.025),
                 ElectricCharge=5,
                 techrequired="unmannedTech")


def _loadAntennas():
    Antenna(name="commDish",
            title="Communotron 88-88",
            mass=0.025,
            price=1100,
            techrequired="electronics")
    Antenna(name="longAntenna",
            title="Communotron 16",
            mass=0.005,
            price=300,
            techrequired="engineering101")
    Antenna(name="mediumDishAntenna",
            title="Comms DTS-M1",
            mass=0.03,
            price=600,
            techrequired="basicScience")


_load()