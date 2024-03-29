from sympy import *
from copy import copy
from time import sleep
from sympy.core.numbers import Float, Integer, Rational, Zero, One, NegativeOne, Half, Exp1, pi, Catalan, EulerGamma, GoldenRatio, TribonacciConstant
from enum import Enum, auto
from pint import UnitRegistry, DimensionalityError

##SIMPLE NOTES
#Paczynskis approx: q is the primary (who we're finding Rl for)/secondary
#Extinction for magnitude is expected - observed

stefanboltzmannconstant = 5.67e-8 #W m-2 K-4
sunmass = 1.99e30 #kg
sunradius = 6.96e5 #km
sunluminosity = 3.9e26 #W
suntemperature = 5778 #K
gasconstant = 8.31446261815324 #Pa m3 K-1 mol-1
gravityconstant = 6.67e-11 #N m2 kg-2
lightspeed = 3e8 #m s-1
hubbleconstant = 70 #km s-1 Mpc-1
hubbletime = 9.778e11/hubbleconstant #yr
earthmass = 5.97e24 #kg
earthradius = 6378 #km
moonmass = 7.346e22 #kg
moonradius = 1738.1 #km
protonmass = 1.6726e-27 #kg
electronmass = 9.109e-31 #kg
planckconstant = 6.626e-34 #J s
boltzmannconstant = 1.380e-23 #m2 kg s-2 K-1
Halphaspectralline = 656.28 #nm
Type_Ia_Supernoval_Absolute_Magnitude = -19.3
sundensity = 1410 #kg m-3
jupiterradius = 69,911 #km
jupitermass = 1.8982e27 #kg

ur = UnitRegistry(system='SI')
ur.autoconvert_offset_to_baseunit=True

ur.define("solar_mass = {} * kg = sunmass = solar_masses = sunmasses".format(sunmass))
ur.define("solar_luminosity = {} * watt = sunluminosity = solar_luminosities = sunluminosities".format(sunluminosity))
ur.define("hubble_time = {} * yr = hubbletime = hubble_times = hubbletimes".format(hubbletime))
ur.define("earth_mass = {} * kg = earthmass = earth_masses = earthmasses".format(earthmass))
ur.define("earth_radius = {} * km = earthradius = earth_radiuses = earthradiuses".format(earthradius))
ur.define("sun_radius = {} * km = sunradius = sun_radiuses = sunradiuses".format(sunradius))
ur.define("sun_temperature = {} * K = suntemperature = sun_temperatures = suntemperatures = sun_temp = suntemp = sun_temps = suntemps".format(suntemperature))
ur.define("moon_mass = {} * kg = moonmass = moon_masses = moonmasses".format(moonmass))
ur.define("moon_radius = {} * km = moonradius = moon_radiuses = moonradiuses".format(moonradius))
ur.define("solar_density = {} * kg * meter**-3 = sundensity = sundensities = sun_densities".format(sundensity))
ur.define("jansky = 1e-26 * watt * meter**-2 * second**-1 = janskies = Jy")
ur.define("jupiter_mass = {} * kg = jupitermass = jupiter_masses = jupitermasses".format(jupitermass))
ur.define("jupiter_radius = {} * km = jupiterradius = jupiter_radiuses = jupiterradiuses".format(jupiterradius))

def pint_val(num):
    return float(num.magnitude) if type(num) not in [int, float, str] else float(num)

class Measure(Enum):
    LENGTH = auto()
    TIME = auto()
    LUMINOSITY = auto()
    FREQUENCY = auto()
    VELOCITY = auto()
    PRESSURE = auto()
    MASS = auto()
    AREA = auto()
    VOLUME = auto()
    QUANTITY = auto()
    FLUX = auto()
    TEMPERATURE = auto()
    ANGLE = auto()
    ENERGY = auto()
    FORCE = auto()
    SQUARE_ANGLE = auto()
    DENSITY = auto()

standard_unit = {
    Measure.LENGTH: ur.meter,
    Measure.TIME: ur.second,
    Measure.LUMINOSITY: ur.watt,
    Measure.FREQUENCY: ur.hertz,
    Measure.VELOCITY: ur.meter_per_second,
    Measure.PRESSURE: ur.atmosphere,
    Measure.MASS: ur.kilogram,
    Measure.AREA: ur.meter**2,
    Measure.VOLUME: ur.meter**3,
    Measure.QUANTITY: ur.mole,
    Measure.FLUX: ur.watt/ur.meter**2,
    Measure.TEMPERATURE: ur.kelvin,
    Measure.ANGLE: ur.radian,
    Measure.ENERGY: ur.joule,
    Measure.FORCE: ur.newton,
    Measure.SQUARE_ANGLE: ur.arcmin**2,
    Measure.DENSITY: ur.kilogram * ur.meter**-3
}
all_units = {i: list(ur.get_compatible_units(standard_unit[i])) for i in standard_unit}
all_units[Measure.MASS] += [ur.solar_mass, ur.earth_mass, ur.moon_mass, ur.proton_mass, ur.electron_mass, ur.jupitermass]
all_units[Measure.LUMINOSITY] += [ur.solar_luminosity]
all_units[Measure.TIME] += [ur.hubble_time]
all_units[Measure.LENGTH] += [ur.earth_radius, ur.sun_radius, ur.moon_radius, ur.nanometer, ur.Mpc, ur.jupiterradius, ur.kilometer]
all_units[Measure.TEMPERATURE] += [ur.sun_temperature]
all_units[Measure.SQUARE_ANGLE] = [ur.steradian, ur.arcsec ** 2, ur.arcmin ** 2, ur.degree ** 2, ur.turn ** 2]
all_units[Measure.DENSITY] += [ur.solar_mass * ur.meter**-3, ur.g * ur.centimeter**-3, ur.g * ur.liter**-1, ur.solar_density, ur.kg * ur.m**-3]
all_units[Measure.VOLUME] += [ur.m**3]

Rp, Ra, a, e = symbols("Rp Ra a e") # aphelion/perihelion
ma, mb, p = symbols("ma mb p") # kepler's third law
d = symbols("d") # simple distance
apmaga, abmaga = symbols("apmaga abmaga") # distance modulus
apmagb, abmagb = symbols("apmagb abmagb") # distance modulus
F, T, L, r = symbols("F T L r") # stefan-boltzmann
lp = symbols("lp") # wien
P, V, n = symbols("P V n") #PVNRT not used yet
Ia, Ib = symbols("Ia Ib") # intensities
thet, D = symbols("thet D") # small angle formula
v = symbols("v") # velocity
lgpa, lgpb, Db = symbols("lgpa lgpb Db") # light gathering power
magpow, Fo, Fe = symbols("magpow Fo Fe") # magnification
lo, le, z = symbols("lo le z") # redshift
E = symbols("E") # E=mc2
pl = symbols("pl") # parallax
fr = symbols("fr") # focal ratio
c = symbols("c") # eccentricity
tau = symbols("tau") # stellar lifetime
Gf = symbols("Gf") # gravitational field
f = symbols("f") # vacuum frequency
C, A, al, h, V, SA = symbols("C A al h V SA") # circumference area perimeter arclength height volume surfacearea
dens = symbols("dens") # density
Fdens = symbols("Fdens") # fluxdensity
densb = symbols("densb") # roche limit
ppr, lpr, epr = symbols("ppr lpr epr")

kep1 = Eq(Rp, a * (1 + e))
kep1 = kep1.subs([(a, 1), (e, 1)])

class EquationRegistry(object):
    def __init__(self, equations):
        self.equations = equations #parallel of unkowns, sympy Eq for each equation

    def append(self, equation, unknown): #unknwown form same as symbols example: "Rp Ra a e"
        self.equations.append(equation)
    
    def solvable(self, knowns): #known form same as symbols
        if type(knowns) is dict:
            knowns = get_knowns(knowns)
        k = knowns.split(" ")
        return [n[1] for n in filter(lambda f: len(f[0]) == 1, [(list(filter(lambda e: e not in k, i[0].split(" "))), i[1]) for i in zip([g.unknowns for g in self.equations], self.equations)])]

class Equation(object):
    def __init__(self, equation, unknowns, name):
        self.name = name
        self.unknowns = unknowns
        self.eq = equation
        self.latex = latex(equation)
    
    def solve(self, values):
        return self.subs_solve([(i, values[i]) for i in self.unknowns.split(" ")])

    def subs_solve(self, subs):
        var = list(filter(lambda e: subs[[i[0] for i in subs].index(e)][1] is None, self.unknowns.split(" ")))[0]
        subs = list(filter(lambda e: e[1] is not None, subs))
        return [n * (units[var] if units[var] is not None else 1) for n in [float(i) for i in list(filter(lambda e: type(e) in [Float, Integer, Rational, Zero, One, NegativeOne, Half, Exp1, pi, Catalan, EulerGamma, GoldenRatio, Integers, TribonacciConstant], solve(self.eq.subs(subs))))]], var
    
    def get_unknowns(self, vals):
        u = []
        for i in self.unknowns.split(" "):
            if vals[i] is None:
                u.append(i)
        return u
    
    def __str__(self):
        return self.name

def is_identity(eq):
    return not any([str(i) in str(eq.lhs) for i in range(10)]) and any([str(i) in str(eq.rhs) for i in range(10)])

def get_knowns(vals):
    out = ""
    for i in vals:
        if vals[i] is not None:
            out += (" " + i) if len(out) > 0 else i
    return out

orbital = [
    Equation(Eq(Ra, a * (1 + e)), "Ra a e", "Aphelion"),
    Equation(Eq(Rp, a * (1 - e)), "Rp a e", "Perihelion"),
    Equation(Eq(ma + mb, (a*6.68459e-12)**3/p**2), "ma mb p a", "Kepler's Third Law"),
    Equation(Eq(v, sqrt(gravityconstant * (ma*sunmass)/r)), "v ma r", "Circular Velocity"),
    Equation(Eq(v, sqrt(gravityconstant * (mb*sunmass)/r)), "v mb r", "Circular Velocity"),
    Equation(Eq(e, c/a), "e c a", "Orbit Eccentricity"),
    Equation(Eq(v, 2 * pi * a / (p*31557600)), "v a p", "Orbital Velocity")
]

stellar = [
    Equation(Eq(L/sunluminosity, ma**3.5), "L ma", "Mass-Luminosity Relation"),
    Equation(Eq(L/sunluminosity, mb**3.5), "L mb", "Mass-Luminosity Relation"),
    Equation(Eq(tau, 10**10 * ma**-2.5), "tau ma", "Stellar Lifetime"),
    Equation(Eq(tau, 10**10 * mb**-2.5), "tau mb", "Stellar Lifetime"),
    Equation(Eq(r, 2 * gravityconstant * (ma*sunmass)/lightspeed**2), "r ma", "Schwarzschild Radius"),
    Equation(Eq(r, 2 * gravityconstant * (mb*sunmass)/lightspeed**2), "r mb", "Schwarzschild Radius")
]

magnitude = [
    Equation(Eq(apmaga - abmaga, -5 + 5 * log(d, 10)), "apmaga abmaga d", "Distance Modulus"),
    Equation(Eq(apmagb - abmagb, -5 + 5 * log(d, 10)), "apmagb abmagb d", "Distance Modulus"),
    Equation(Eq(apmagb - apmaga, 2.5 * log(Ia/Ib, 10)), "apmagb apmaga Ia Ib", "Intensity Ratio"),
    Equation(Eq(Ia, L/(4 * pi * (d*3.086e16)**2)), "Ia L d", "Inverse Squares Law"),
    Equation(Eq(Ib, L/(4 * pi * (d*3.086e16)**2)), "Ib L d", "Inverse Squares Law"),
    Equation(Eq(abmaga, -2.5 * log(L*3.0128e-28, 10)), "abmaga L", "Magnitude Luminosity Relationship"),
    Equation(Eq(abmagb, -2.5 * log(L*3.0128e-28, 10)), "abmagb L", "Magnitude Luminosity Relationship")
]

stefan_boltzmann = [
    Equation(Eq(F, stefanboltzmannconstant * T**4), "F T", "Stefan-Boltzmann Law"),
    Equation(Eq(L, 4 * pi * stefanboltzmannconstant * r**2 * T**4), "L r T", "Stefan-Boltzmann Law"),
    Equation(Eq(L/sunluminosity, (r/sunradius)**2 * (T/suntemperature)**4), "L r T", "Stefan-Boltzmann Law"),
    Equation(Eq(Fdens, F * f), "F f Fdens", "Flux Density")
]

wavelength = [
    Equation(Eq(lp, 2900000/T), "lp T", "Wien's Law"),
    Equation(Eq(z + 1, lo/le), "z lo le", "Redshift"),
    Equation(Eq(v, lightspeed * z), "v z", "Doppler Formula"),
    Equation(Eq(z + 1, sqrt((1 + (v/lightspeed))/(1 + (v/lightspeed)))), "v z", "Relativistic Doppler Shift"),
    Equation(Eq(f, lightspeed/(lp/1e9)), "f lp", "Vacuum Frequency"),
    Equation(Eq(f, lightspeed/(lo/1e9)), "f lo", "Vacuum Frequency"),
    Equation(Eq(f, lightspeed/(le/1e9)), "f le", "Vacuum Frequency"),
]

gas = [
    Equation(Eq(P * V, n * gasconstant * T), "P V n T", "Ideal Gas Law"),
]

relativity = [
    Equation(Eq(E, ma/sunmass * lightspeed**2), "E ma", "Fusion"),
    Equation(Eq(E, mb/sunmass * lightspeed**2), "E mb", "Fusion"),
    Equation(Eq(ppr, p/sqrt(1 - (v**2/lightspeed**2))), "ppr p v", "Time Dilation"),
    Equation(Eq(lpr, D * sqrt(1 - (v**2/lightspeed**2))), "lpr D v", "Length Contraction"),
    Equation(Eq(lpr, h * sqrt(1 - (v**2/lightspeed**2))), "lpr h v", "Length Contraction"),
    Equation(Eq(epr, (ma * lightspeed**2)/sqrt(1 - (v**2/lightspeed**2))), "epr ma v", "Relativistic Energy"),
    Equation(Eq(epr, (mb * lightspeed**2)/sqrt(1 - (v**2/lightspeed**2))), "epr mb v", "Relativistic Energy")
]

angles_and_telescopes = [
    Equation(Eq(thet, 206265 * D/d), "thet D d", "Small Angle Formula"),
    Equation(Eq(thet, 251643 * (lp/1e9)/D), "thet lp D", "Resolving Power"), #lp is just lambda here, D is optic diameter
    Equation(Eq(thet, 251643 * (lo/1e9)/D), "thet lo D", "Resolving Power"), #lo is just lambda here,
    Equation(Eq(thet, 251643 * (le/1e9)/D), "thet le D", "Resolving Power"), #lp is just lambda here,
    Equation(Eq(lgpa/lgpb, (D/Db)**2), "lgpa lgpb D Db", "Compare LGP"),
    Equation(Eq(magpow, Fo/Fe), "magpow Fo Fe", "Magnification"),
    Equation(Eq(d, 1/pl), "d pl", "Parallax"),
    Equation(Eq(fr, Fe/D), "fr Fe D", "Focal Ratio"),
    Equation(Eq(fr, Fo/D), "fr Fo D", "Focal Ratio")
]

recession = [
    Equation(Eq(v/1000, hubbleconstant * (d*1e-6)), "v d", "Hubble's Law")
]

gravitation = [
    Equation(Eq(Gf, gravityconstant * ma * mb/r**2), "Gf ma mb r", "Universal Gravitation"),
    Equation(Eq((d*3.08567758e+13), 2.4 * (r/1000) * (dens/densb)**(1/3)), "dens densb r d", "Roche Limit")
]

geometry = [
    Equation(Eq(D, 2 * r), "D r", "Diameter Radius Relation"),
    Equation(Eq(C, 2 * pi * r), "C r", "Circumference of a Circle"),
    Equation(Eq(A, pi * r**2), "A r", "Area of a Circle"),
    Equation(Eq((thet*pi/(180*3600))/(2 * pi), al/C), "thet al C", "Arc Length"),
    Equation(Eq(A, (thet*pi/(180*3600)) * r**2 * .5), "thet r A", "Arc Area"),
    Equation(Eq(A * h, V), "A h V", "Volume of a Cylinder"),
    Equation(Eq(SA, 2 * A + C * h), "SA A C h", "Surface Area of a Cylinder"),
    Equation(Eq(SA, pi * r * (r + sqrt(h**2 + r**2))), "SA r h", "Surface Area of a Cone"),
    Equation(Eq(V, 1/3 * pi * r**2 * h), "V r h", "Volume of a Cone"),
    Equation(Eq(V, 4/3 * pi * r**3), "V r", "Volume of a Sphere"),
    Equation(Eq(SA, 4 * pi * r**2), "SA r", "Surface Area of a Sphere"),
    Equation(Eq(dens, ma/V), "dens ma V", "Density"),
    Equation(Eq(dens, mb/V), "dens mb V", "Density")
]

all_equations = orbital + magnitude + stefan_boltzmann + wavelength + gas + angles_and_telescopes + recession + stellar + gravitation + geometry + relativity

er = EquationRegistry(all_equations)

vals = {
    "Rp" : None,
    "Ra" : None,
    "a" : None,
    "e" : None,
    "apmaga": None,
    "abmaga": None,
    "apmagb": None,
    "abmagb": None,
    "d": None,
    "p": None,
    "ma": None,
    "mb": None,
    "F": None,
    "T": None,
    "L": None,
    "r": None,
    "P": None,
    "V": None,
    "n": None,
    "lp": None,
    "Ia": None,
    "Ib": None,
    "thet": None,
    "D": None,
    "v": None,
    "lgpa": None,
    "lgpb": None,
    "Db": None,
    "magpow": None,
    "Fo": None,
    "Fe": None,
    "lo": None,
    "le": None,
    "z": None,
    "E": None,
    "pl": None, 
    "fr": None,
    "c": None,
    "tau": None,
    "Gf": None,
    "f": None,
    "C": None,
    "A": None,
    "al": None,
    "h": None,
    "SA": None,
    "dens": None,
    "Fdens": None,
    "densb": None,
    "ppr": None,
    "lpr": None,
    "epr": None
}

units = {
    "Rp" : ur.meters,
    "Ra" : ur.meters,
    "a" : ur.meters,
    "e" : None,
    "apmaga": None,
    "abmaga": None,
    "d": ur.parsecs,
    "apmagb": None,
    "abmagb": None,
    "ma": ur.solar_mass,
    "mb": ur.solar_mass,
    "p": ur.year,
    "F": ur.watts * ur.meters**-2,
    "T": ur.kelvins,
    "L": ur.watts,
    "r": ur.meters,
    "P": ur.pascals,
    "V": ur.meters**3,
    "n": ur.moles,
    "lp": ur.nanometers,
    "Ia": ur.watts * ur.meters**-2,
    "Ib": ur.watts * ur.meters**-2,
    "thet": ur.arcsec,
    "D": ur.meters,
    "v": ur.meters * ur.second**-1,
    "lgpa": None,
    "lgpb": None,
    "Db": ur.meters,
    "magpow": None,
    "Fo": ur.meters,
    "Fe": ur.meters,
    "lo": ur.nanometers,
    "le": ur.nanometers,
    "z": None,
    "E": ur.joules,
    "pl": ur.arcsec, 
    "fr": None,
    "c": ur.meters,
    "tau": ur.years,
    "Gf": ur.newtons,
    "f": ur.second**-1,
    "C": ur.meters,
    "A": ur.meters**2,
    "al": ur.meters,
    "h": ur.meters,
    "SA": ur.meters**2,
    "dens": ur.solar_mass * ur.meters**-3,
    "Fdens": ur.jansky,
    "densb": ur.solar_mass * ur.meters**-3,
    "ppr": ur.seconds,
    "lpr": ur.meters,
    "epr": ur.joules
}

process = {}
def retlamb1():
    return lambda e: e
def retlamb2(i):
    return lambda e: e if type(e) in [int, float] else pint_val(ur.parse_expression(e).to(units[i]))
for i in units:
    if units[i] is None:
        process[i] = retlamb1()
    else:
        process[i] = retlamb2(i)

def set_val(key, val):
    if val is not None:
        vals[key] = process[key](val)
    else:
        vals[key] = None

def pre_screen(key, val):
    try:
        process[key](val)
    except:
        return False
    return True

def as_all_units(var, num, meas="nomeas"):
    if meas == "nomeas":
        meas = measure_of[var]
    if meas is not None:
        ret = [num]
        for i in all_units[meas]:
            try:
                ret.append(num.to(i))
            except Exception as err:
                print(type(err), err)
        return ret
    else:
        return [num]

long_name = {
    "Rp" : "Orbit Perihelion",
    "Ra" : "Orbit Aphelion",
    "a" : "Semimajor Axis",
    "e" : "Eccentricity",
    "apmaga": "Apparent Magnitude of Object A",
    "abmaga": "Absolute Magnitude of Object A",
    "d": "Distance",
    "apmagb": "Apparent Magnitude of Object B",
    "abmagb": "Absolute Magnitude of Object B",
    "ma": "Mass of Object A",
    "mb": "Mass of Object B",
    "p": "Orbital Period",
    "F": "Flux",
    "T": "Temperature",
    "L": "Luminosity",
    "r": "Radius",
    "P": "Pressure",
    "V": "Volume",
    "n": "Quantity",
    "lp": "Peak Wavelength",
    "Ia": "Intensity of Object A",
    "Ib": "Intensity of Object B",
    "thet": "Theta",
    "D": "Diameter",
    "v": "Velocity",
    "lgpa": "Light Gathering Power of Object A",
    "lgpb": "Light Gathering Power of Object B",
    "Db": "Diameter of Object B",
    "magpow": "Magnification Power",
    "Fo": "Objective Focal Length",
    "Fe": "Eyepiece Focal Length",
    "lo": "Observed Wavelength",
    "le": "Emitted Wavelength",
    "z": "Redshift",
    "E": "Energy",
    "pl": "Parallax", 
    "fr": "Focal Ratio",
    "c": "Distance from Focus to Center",
    "tau": "Stellar Lifetime",
    "Gf": "Gravitational Force",
    "C": "Circumeference",
    "A": "Area",
    "al": "Arc Length",
    "h": "Height",
    "SA": "Surface Area",
    "dens": "Density",
    "f": "Frequency",
    "Fdens": "Flux Density",
    "densb": "Density of Object B",
    "ppr": "Contracted Time",
    "lpr": "Contracted Length",
    "epr": "Relative Energy"
}

measure_of = {
    "Rp": Measure.LENGTH,
    "Ra": Measure.LENGTH,
    "a": Measure.LENGTH,
    "e": None,
    "apmaga": None,
    "abmaga": None,
    "apmagb": None,
    "abmagb": None,
    "ma": Measure.MASS,
    "mb": Measure.MASS,
    "d": Measure.LENGTH,
    "p": Measure.TIME,
    "F": Measure.FLUX,
    "T": Measure.TEMPERATURE,
    "L": Measure.LUMINOSITY,
    "r": Measure.LENGTH,
    "P": Measure.PRESSURE,
    "V": Measure.VOLUME,
    "n": Measure.QUANTITY,
    "lp": Measure.LENGTH,
    "Ia": Measure.FLUX,
    "Ib": Measure.FLUX,
    "thet": Measure.ANGLE,
    "D": Measure.LENGTH,
    "v": Measure.VELOCITY,
    "lgpa": None,
    "lgpb": None,
    "Db": Measure.LENGTH,
    "magpow": None,
    "Fo": Measure.LENGTH,
    "Fe": Measure.LENGTH,
    "lo": Measure.LENGTH,
    "le": Measure.LENGTH,
    "z": None,
    "E": Measure.ENERGY,
    "pl": Measure.ANGLE, 
    "fr": None,
    "c": Measure.LENGTH,
    "tau": Measure.TIME,
    "Gf": Measure.FORCE,
    "C": Measure.LENGTH,
    "A": Measure.AREA,
    "al": Measure.LENGTH,
    "h": Measure.LENGTH,
    "SA": Measure.AREA,
    "dens": Measure.DENSITY,
    "f": Measure.FREQUENCY,
    "Fdens": None,
    "densb": Measure.DENSITY,
    "ppr": Measure.TIME,
    "lpr": Measure.LENGTH,
    "epr": Measure.ENERGY
}

if __name__ == "__main__":
    print("Rp, a", er.solvable("Rp a"))
    print("Rp, a, e", er.solvable("Rp a e"))
    print("Rp, Ra, e", er.solvable("Rp Ra e"))

    vals["Rp"] = 14637
    vals["a"] = 2
    print(er.solvable("Rp a")[0].subs_solve([(Rp, vals["Rp"]), (a, vals["a"])]))
    print(er.solvable("Rp a")[0].solve(vals))