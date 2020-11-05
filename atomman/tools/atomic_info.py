# coding: utf-8

__all__ = ['atomic_number', 'atomic_symbol', 'atomic_mass']

import numpy as np
import pandas as pd

from .atomicdata import atomicdata

class AtomicInfo():
    
    def __init__(self):
        """Class initializer"""
        self.__data = atomicdata
    
    @property
    def data(self):
        return self.__data
    
    @property
    def renames(self):
        """dict : Matches systematic element symbols to their now assigned symbols"""
        return {
            'Unq': 'Rf',
            'Unp': 'Db',
            'Unh': 'Sg',
            'Uns': 'Bh',
            'Uno': 'Hs',
            'Une': 'Mt',
            'Uun': 'Ds',
            'Uuu': 'Rg',
            'Uub': 'Cn',
            'Uut': 'Nh',
            'Uuq': 'Fl',
            'Uup': 'Mc',
            'Uuh': 'Lv',
            'Uus': 'Ts',
            'Uuo': 'Og',
        }        
    
    def load(self, datafile):
        """
        Reads in the "Linearized ASCII Output" as found at https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
        """
        with open(datafile) as f:
            lines = f.readlines()

        data = []
        i = 0
        while i < len(lines):
            if len(lines[i]) > 15 and lines[i][:15] == 'Atomic Number =':
                d = {}
                while i < len(lines) and len(lines[i].strip()) > 0:
                    terms = lines[i].split('=')
                    if len(terms) > 1:
                        d[terms[0].strip()] = terms[1].strip()
                    else:
                        d[terms[0].strip()] = np.nan 
                    i += 1
                data.append(d)
            else:
                i += 1
        data = pd.DataFrame(data)
        
        def make_int(series, key):
            return int(series[key])
        
        data['Atomic Number'] = data.apply(make_int, args=['Atomic Number'], axis=1)
        data['Mass Number'] = data.apply(make_int, args=['Mass Number'], axis=1)
        
        self.__data = data
        
    def atomic_number(self, atomic_symbol):
        """
        Return the corresponding atomic number for a given atomic symbol.

        Parameters
        ----------
        atomic_symbol : str
            An atomic symbol.

        Return
        ------
        int
            The corresponding atomic number.
        """
        
        # Handle old systematic named symbols
        if atomic_symbol in self.renames:
            atomic_symbol = self.renames[atomic_symbol]
        
        matches = self.data[self.data['Atomic Symbol'] == atomic_symbol]
        if len(matches) > 0:
            return matches.iloc[0]['Atomic Number']
        else:
            raise ValueError(f'No matches for atomic symbol {atomic_symbol} found')
    
    def atomic_symbol(self, atomic_number):
        """
        Return the corresponding atomic symbol for a given atomic number.

        Parameters
        ----------
        atomic_number : int
            An atomic number.

        Return
        ------
        int
            The corresponding atomic symbol.
        """
        matches = self.data[self.data['Atomic Number'] == atomic_number]
        if len(matches) > 0:
            return matches.iloc[0]['Atomic Symbol']
        else:
            raise ValueError(f'No matches for atomic number {atomic_number} found')
    
    def atomic_mass(self, atomic_info, mass_number=None, prompt=True):
        """
        Returns either the median standard atomic weight for an element or the relative
        atomic mass for an isotope.
        
        Parameters
        ----------
        atomic_info : str or int
            The atomic symbol or number identifying the element/isotopes.  
        mass_number : int, optional
            An isotope mass number.
        Returns
        -------
        float
            The average standard atomic weight of an element or the relative
            atomic mass of an isotope.
        """

        # Try converting atomic_info to an int - fetch atomic_symbol if needed
        try:
            atomic_number = int(atomic_info)
        except:
            atomic_symbol = atomic_info
        else:
            atomic_symbol = self.atomic_symbol(atomic_number)        
            
        # Separate mass_number from atomic_symbol if needed
        if '-' in atomic_symbol:
            if mass_number is not None:
                raise ValueError('Mass number cannot be given both with the symbol and separately')
            atomic_symbol, mass_number = atomic_symbol.split('-')
        
        # Handle hydrogens
        if atomic_symbol in ['H', 'D', 'T']:
            atomic_symbol, mass_number = self.__handle_hydrogen(atomic_symbol, mass_number=mass_number)
        
        # Handle old systematic named symbols
        if atomic_symbol in self.renames:
            atomic_symbol = self.renames[atomic_symbol]
        
        # Check if there is a standard atomic weight for an element
        if mass_number is None:
            matches = self.data[self.data['Atomic Symbol'] == atomic_symbol]
            weight = matches.iloc[0]['Standard Atomic Weight']
            if isinstance(weight, float):
                if pd.notna(weight):
                    return weight
            elif '[' in weight:
                return np.array(eval(weight)).mean()
            elif '(' in weight:
                return float(weight.split('(')[0])
            elif weight != '':
                return float(weight)
            elif len(matches) == 0:
                raise ValueError(f'No matches for atomic symbol {atomic_symbol} found')
        
            # Return isotope mass if only one isotope
            if len(matches) == 1:
                return float(matches.iloc[0]['Relative Atomic Mass'].split('(')[0])
            
            if prompt:
                print(f'No standard atomic weight for {atomic_symbol}.')
                isotopes = matches['Mass Number'].to_list()
                print(f'Please select an isotope from {isotopes}:')
                mass_number = input()
            else:
                raise ValueError(f'No standard weight for {atomic_symbol}. Specify an isotope instead.')
        
        # Convert mass_number to an int
        mass_number = int(mass_number)
                      
        # Find relative atomic mass
        matches = self.data[(self.data['Atomic Symbol'] == atomic_symbol) 
                           &(self.data['Mass Number'] == mass_number)]
        
        if len(matches) == 1:
            mass = matches.iloc[0]['Relative Atomic Mass']
            if isinstance(mass, float):
                return mass
            elif '(' in mass:
                return float(mass.split('(')[0])
            else:
                raise ValueError('Mass value format not recognized!!!!!!!!!!')

        elif len(matches) == 0:
            raise ValueError(f'No matches for atomic symbol {atomic_symbol} and mass number {mass_number} found')
        else:
            raise ValueError('Multiple matches found!!!')
        
    def __handle_hydrogen(self, atomic_symbol, mass_number=None):
        """
        Special handling of hydrogen isotopes due to isotope symbols D and T
        
        Parameters
        ----------
        atomic_symbol : str
            The atomic symbol.  
        mass_number : int, optional
            An isotope mass number.
        
        Returns
        -------
        atomic_symbol : str
            The correct database atomic symbol for the input
        mass_number : int or None
            The correct database mass number for the input
        """
        if mass_number is not None:
            mass_number = int(mass_number)
            
        # Check that mass number is 2 for 'D'
        if atomic_symbol == 'D':
            if mass_number is None:
                mass_number = 2
            elif mass_number != 2:
                raise ValueError('Mass number must be 2 for isotope symbol D')
        
        # Check that mass number is 3 for 'T'
        elif atomic_symbol == 'T':
            if mass_number is None:
                mass_number = 3
            elif mass_number != 3:
                raise ValueError('Mass number must be 3 for isotope symbol T')
        else:
            # Change 'H' to 'D' for mass number 2
            if mass_number == 2:
                atomic_symbol = 'D'
            # Change 'H' to 'T' for mass number 3
            elif mass_number == 3:
                atomic_symbol = 'T'
                
        return atomic_symbol, mass_number

# Initialize an object of the class and create shortcuts to the methods
atomicinfo = AtomicInfo()
atomic_number = atomicinfo.atomic_number
atomic_symbol = atomicinfo.atomic_symbol
atomic_mass = atomicinfo.atomic_mass