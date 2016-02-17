#External library imports
import numericalunits as nu
import numpy as np

def build_unit():
    """Builds the unit dictionary based on current working units"""
    
    global unit
    unit = {
        "m": nu.m, "cm": nu.cm, "mm": nu.mm, "um": nu.um, "nm": nu.nm, "pm": nu.pm, "fm": nu.fm, 
        "km": nu.km, "angstrom": nu.angstrom, "Angstrom": nu.angstrom, "lightyear": nu.lightyear, "astro_unit": nu.astro_unit,
        "pc": nu.pc, "kpc": nu.kpc, "Mpc": nu.Mpc, "Gpc": nu.Gpc, 
        "inch": nu.inch, "foot": nu.foot, "mile": nu.mile, "thou": nu.thou,
        
        "L": nu.L, "mL": nu.mL, "uL": nu.uL, "nL": nu.nL, "pL": nu.pL, "fL": nu.fL,
        "aL": nu.aL, "kL": nu.kL, "ML": nu.ML, "GL": nu.GL,

        "s": nu.s, "ms": nu.ms, "us": nu.us, "ns": nu.ns, "ps": nu.ps, "fs": nu.fs,
        "minute": nu.minute, "hour": nu.hour, "day": nu.day, "week": nu.week, "year": nu.year,
        
        "Hz": nu.Hz, "mHz": nu.mHz, "kHz": nu.kHz, "MHz": nu.MHz, "GHz": nu.GHz,
        "THz": nu.THz, "PHz": nu.PHz,
        
        "kg": nu.kg, "g": nu.g, "mg": nu.mg, "ug": nu.ug, "ng": nu.ng, "pg": nu.pg, "fg": nu.fg,
        "tonne": nu.tonne, "amu": nu.amu,
        "Da": nu.Da, "kDa": nu.kDa, "lbm": nu.lbm, 
        
        "J": nu.J, "mJ": nu.mJ, "uJ": nu.uJ, "nJ": nu.nJ, "pJ": nu.pJ, "fJ": nu.fJ, 
        "kJ": nu.kJ, "MJ": nu.MJ, "GJ": nu.GJ, "erg": nu.erg,
        "eV": nu.eV, "meV": nu.meV, "keV": nu.keV, "MeV": nu.MeV, "GeV": nu.GeV,
        "TeV": nu.TeV, "btu": nu.btu, "smallcal": nu.smallcal, "kcal": nu.kcal,
        "Wh": nu.Wh, "kWh": nu.kWh,
        
        "NA": nu.NA, "mol": nu.mol, "mmol": nu.mmol, "umol": nu.umol, "nmol": nu.nmol,
        "pmol": nu.pmol, "fmol": nu.fmol,
        "M": nu.M, "mM": nu.mM, "uM": nu.uM, "nM": nu.nM, "pM": nu.pM, "fM": nu.fM,
        
        "N": nu.N, "dyn": nu.dyn, "lbf": nu.lbf,

        "Pa": nu.Pa, "hPa": nu.hPa, "kPa": nu.kPa, "MPa": nu.MPa, "GPa": nu.GPa,  
        "bar": nu.bar, "mbar": nu.mbar, "cbar": nu.cbar, "dbar": nu.dbar, "kbar": nu.kbar, "Mbar": nu.Mbar,
        "atm": nu.atm, "torr": nu.torr, "mtorr": nu.mtorr, "psi": nu.psi,
        
        "W": nu.W, "mW": nu.mW, "uW": nu.uW, "nW": nu.nW, "pW": nu.pW, "kW": nu.kW,
        "MW": nu.MW, "GW": nu.GW, "TW": nu.TW,

        "K": nu.K, "degFinterval": nu.degFinterval, "degCinterval": nu.degCinterval,
        
        "C": nu.C, "mC": nu.mC, "uC": nu.uC, "nC": nu.nC, "Ah": nu.Ah, "mAh": nu.mAh,
        
        "A": nu.A, "mA": nu.mA, "uA": nu.uA, "nA": nu.nA, "pA": nu.pA, "fA": nu.fA,
        
        "V": nu.V, "mV": nu.mV, "uV": nu.uV, "nV": nu.nV, "kV": nu.kV, "MV": nu.MV, "GV": nu.GV, "TV": nu.TV,
        
        "ohm": nu.ohm, "mohm": nu.mohm, "kohm": nu.kohm, "Mohm": nu.Mohm, "Gohm": nu.Gohm,
        
        "S": nu.S, "mS": nu.mS, "uS": nu.uS, "nS": nu.nS,
        
        "T": nu.T, "mT": nu.mT, "uT": nu.uT, "nT": nu.nT,
        "G": nu.G, "mG": nu.mG, "uG": nu.uG, "kG": nu.kG,
        "Oe": nu.Oe, "Wb": nu.Wb,
        
        "F": nu.F, "uF": nu.uF, "nF": nu.nF, "pF": nu.pF, "fF": nu.fF, "aF": nu.aF,
        
        "H": nu.H, "mH": nu.mH, "uH": nu.uH, "nH": nu.nH,

        "c0": nu.c0, "mu0": nu.mu0, "eps0": nu.eps0, "Z0": nu.Z0, "hPlanck": nu.hPlanck, "hbar": nu.hbar,
        "kB": nu.kB, "GNewton": nu.GNewton, "sigmaSB": nu.sigmaSB, "alphaFS": nu.alphaFS, "Rgas": nu.Rgas, 
        "e": nu.e, "uBohr": nu.uBohr, "uNuc": nu.uNuc, "aBohr": nu.aBohr, "me": nu.me, "mp": nu.mp, 
        "mn": nu.mn, "Rinf": nu.Rinf, "Ry": nu.Ry, "ARichardson": nu.ARichardson, "Phi0": nu.Phi0,
        "KJos": nu.KJos, "RKlitz": nu.RKlitz, "REarth": nu.REarth, "g0": nu.g0, "Msolar": nu.Msolar,
        "MEarth": nu.MEarth
    }

def reset_units(seed=None, length=None, mass=None, time=None, energy=None, charge=None):
    """Allows the working units to be reset to random values or specified accoring to dimensions."""
    
    if length is None and time is None and mass is None and energy is None and charge is None:
        nu.reset_units(seed)
    
    else:
        try:
            length = unit[length]
        except:
            pass
        try:
            mass = unit[mass]
        except:
            pass
        try:
            time = unit[time]
        except:
            pass        
        try:
            energy = unit[energy]
        except:
            pass    
        try:
            charge = unit[charge]
        except:
            pass      
    
        m = 1.
        kg = 1.
        s = 1.
        C = 1.
        K = 1.
        
        if length is not None:
            m = nu.m / length
        if mass is not None:
            kg = nu.kg / mass
        if time is not None:
            s = nu.s / time
        if charge is not None:
            C = nu.C / charge
        if energy is not None:
            J = nu.J / energy
            if mass is None:
                kg = J * s**2 / m**2
            elif time is None:
                s = (kg * m**2 / J)**0.5
            elif length is None:
                m = (J * s**2 / kg)
            else:
                raise ValueError('length, mass, time and energy cannot all be defined')
        
        nu.m = m
        nu.kg = kg
        nu.s = s
        nu.C = C
        nu.K = K
        
        nu.set_derived_units_and_constants()
    build_unit()    

def set_in_units(value, units):
    """Convert value from specified units to working units"""
    units = parse(units)
    return np.asarray(value) * units
        
def get_in_units(value, units):
    """Convert value from working units to specified units"""
    units = parse(units)
    return np.asarray(value) / units
    
def value_unit(term):
    """Used for data model elements with value, unit"""
    try:
        unit = term['unit']
    except:
        unit = None
    return set_in_units(term['value'], unit)
    
def parse(units):
    """Convert units as strings (or None) into numbers."""
    
    #Units of None does no scaling
    if units is None:
        return 1

    #Parse string and return number value
    elif isinstance(units, (str, unicode)):
        i = 0
        terms = []
        while i < len(units):
            if units[i] == '(':
                end = units.rindex(')')
                terms.append(parse(units[i+1:end]))
                i = end+1
            elif units[i].isalpha():
                term = ''
                while i < len(units) and units[i] not in ' */^\n\r\t':
                    term += units[i]
                    i += 1
                terms.append(unit[term])
            elif units[i].isdigit() or units[i] == '-' or units[i] == '.':
                term = ''
                while i < len(units) and units[i] not in ' */^\n\r\t':
                    term += units[i]
                    i += 1
                terms.append(float(term))
            elif units[i] in '*/^':
                terms.append(units[i])
                i += 1
            elif units[i] in ' \n\r\t':
                i += 1
            else:
                raise ValueError('Unknown character: %s' % units[i])
        
        while '^' in terms:
            c = terms.index('^')
            value = [terms[c-1] ** terms[c+1]]
            terms = terms[:c-1] + value + terms[c+2:]
        
        while len(terms) > 1:
            if terms[1] == '*':
                value = [terms[0] * terms[2]]
                terms = value + terms[3:]
            elif terms[1] == '/':
                value = [terms[0] / terms[2]]
                terms = value + terms[3:]
            else:
                print terms
                raise ValueError('Invalid string format')
        
        return terms[0]

    #Else assume units is already a number
    else:
        return units
   
        
#Initial build and reset.  Only called first time module is loaded        
build_unit()
reset_units()