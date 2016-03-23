from DataModelDict import DataModelDict
from atomman.tools import atomic_mass
import os

class Potential(object):
    """class for building LAMMPS input lines from a LAMMPS-potential data model.""" 
    
    def __init__(self, model, pot_dir=''):
        """
        initializes an instance associated with a LAMMPS-potential data model.
        
        Arguments:
        model -- a string or file-like obect of a json/xml data model containing a LAMMPS-potential branch.
        pot_dir -- (optional) the directory location of any artifacts associated with the potential.
        """
        
        self.load(model, pot_dir)
            
    def load(self, model, pot_dir=None):
        """
        loads LAMMPS-potential data model info.
        
        Arguments:
        model -- a string or file-like obect of a json/xml data model containing a LAMMPS-potential branch.
        pot_dir -- (optional) the directory location of any artifacts associated with the potential.
        """

        self.__dm = DataModelDict(model).find('LAMMPS-potential')        
        
        for atom in self.__dm.iteraslist('atom'):
            #Check if element is listed
            try:
                test = atom['element']
            #If no element is listed, symbol and mass must be
            except:
                try:
                    test = atom['symbol']
                    test = atom['mass']
                    atom['element'] = atom['symbol']
                except:
                    raise KeyError("Error reading Potential's atomic info: mass and symbol are needed if element is not given!")
            
            #Check if symbol is listed.  If not, make symbol = element
            try:
                test = atom['symbol']
            except:
                atom['symbol'] = atom['element']
            
            #Check if mass is listed.  If not, set to standard value of element
            try:
                mass_check = atom['mass']
            except:
                atom['mass'] = atomic_mass(atom['element'])

            assert isinstance(atom['mass'], float), 'Mass needs to be a number!'
        
        if pot_dir is not None:
            self.pot_dir = pot_dir
            
    def __str__(self):
        """The string of the Potential returns its human-readable id"""
        return self.id
    
    @property
    def pot_dir(self):
        """The directory containing files associated with a given potential."""
        return str(self.__pot_dir)
        
    @pot_dir.setter
    def pot_dir(self, value):
        self.__pot_dir = str(value)
        
    @property
    def id(self):
        """Human-readable identifier."""
        return self.__dm['potential']['id']
    
    @property    
    def uuid(self):
        """uuid hash-key."""
        return self.__dm['potential']['key']
    
    @property
    def units(self):
        """LAMMPS units option."""
        return self.__dm['units']
    
    @property
    def atom_style(self):
        """LAMMPS atom_style option."""
        return self.__dm['atom_style']
    
    @property
    def symbols(self):
        """List of all atom-model symbols."""
        symbols = []
        for atom in self.__dm.iteraslist('atom'):
            symbols.append(str(atom['symbol']))    
        return symbols
    
    def elements(self, symbols=None):
        """
        Return list of element names associated with a list of atom-model symbols.
        
        Arguments:
        symbols -- List of atom-model symbols.  If None (default), will use all of the Potential's symbols, i.e. Potential.symbols.
        """
        if symbols is None:
            symbols = self.symbols
        
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]
        
        elements = []
        for symbol in symbols:
            for atom in self.__dm.iteraslist('atom'):
                if symbol == atom['symbol']:
                    elements.append(str(atom['element']))
                    break
        assert len(symbols) == len(elements), 'Not all elements found!'
        
        return elements
            
    def masses(self, symbols=None):
        """
        Return list of element masses associated with a list of atom-model symbols.
        
        Arguments:
        symbols -- List of atom-model symbols.  If None (default), will use all of the Potential's symbols, i.e. Potential.symbols.
        """
        if symbols is None:
            symbols = self.symbols
        
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]
        
        masses = []
        for symbol in symbols:
            for atom in self.__dm.iteraslist('atom'):
                if symbol == atom['symbol']:
                    masses.append(atom['mass'])
                    break
        assert len(symbols) == len(masses), 'Not all masses found!'
        
        return masses
           
    def pair_info(self, symbols = None):
        """
        Return all LAMMPS input command lines associated with the Potential and a list of atom-model symbols.
        
        Arguments:
        symbols -- (optional) list of atom-model symbols being used.  If None (default), will use all of the Potential's elements.
        """
        #if no symbols supplied use all for potential
        if symbols is None:
            symbols = self.symbols
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]
        
        #Generate mass lines
        masses = self.masses(symbols)
        mass = ''
        for i in xrange(len(masses)):
            mass += 'mass %i %f' % ( i+1, masses[i] ) + '\n'
        mass +='\n'
        
        #Generate pair_style line
        style = 'pair_style ' + self.__dm['pair_style']['type'] 
        terms = self.__dm['pair_style'].get('term', None)
        style += self.__pair_terms(self.__dm['pair_style'].iteraslist('term')) + '\n'
       
        #Generate pair_coeff lines
        coeff = ''
        for coeff_line in self.__dm.iteraslist('pair_coeff'):
            if 'interaction' in coeff_line:
                interaction = coeff_line['interaction'].get('symbol', ['*', '*'])            
            else:
                interaction = ['*', '*']
            
            #Always include coeff lines that act on all atoms in the system
            if interaction == ['*', '*']:
                coeff_symbols = self.symbols
                coeff += 'pair_coeff * *' + self.__pair_terms(coeff_line.iteraslist('term'), symbols, coeff_symbols) + '\n'
                continue
            
            #Many-body potentials will contain a symbols term
            if len(coeff_line.finds('symbols')) > 0:
                many = True
            else:
                many = False
                                                    
            #Treat as a many-body potential
            if many:
                coeff_symbols = interaction
                coeff += 'pair_coeff * *' + self.__pair_terms(coeff_line.iteraslist('term'), symbols, coeff_symbols) + '\n'
                
            #Treat as pair potential
            else:
                coeff_symbols = interaction
                assert len(coeff_symbols) == 2,     'Pair potential interactions need two listed elements'
                
                #Classic eam style is a special case
                if self.__dm['pair_style']['type'] == 'eam':
                    assert coeff_symbols[0] == coeff_symbols[1], 'Only i==j interactions allowed for eam style'
                    for i in xrange( len(symbols) ):
                        if symbols[i] == coeff_symbols[0]:
                            coeff += 'pair_coeff %i %i' % (i + 1, i + 1) + self.__pair_terms(coeff_line.iteraslist('term'), symbols, coeff_symbols) + '\n'
                
                #All other pair potentials
                else:
                    for i in xrange( len(symbols) ):
                        for j in xrange( i, len(symbols) ):
                            if (symbols[i] == coeff_symbols[0] and symbols[j] == coeff_symbols[1]) or (symbols[i] == coeff_symbols[1] and symbols[j] == coeff_symbols[0]):
                               
                                coeff += 'pair_coeff %i %i' % (i + 1, j + 1) + self.__pair_terms(coeff_line.iteraslist('term'), symbols, coeff_symbols) + '\n'

        #generate additional command lines
        command = ''
        
        for command_line in self.__dm.iteraslist('command'):
            command += self.__pair_terms(command_line.iteraslist('term'), symbols, self.symbols).strip() + '\n'

        return mass + style + coeff + command
    
    def __pair_terms(self, terms, system_symbols = [], coeff_symbols = []):
        """utility function used by self.pair_info() for composing lines from terms"""
        line = ''
        
        for term in terms:
            for ttype, tval in term.iteritems():
                #print options and parameters as strings
                if ttype == 'option' or ttype == 'parameter': 
                    line += ' ' + str(tval)
                    
                #print files with pot_dir prefixes    
                elif ttype == 'file':
                    line += ' ' + str( os.path.join(self.pot_dir, tval) )

                #print all symbols being used for symbolsList
                elif ttype == 'symbolsList' and (tval is True or tval == 'True'):
                    for coeff_symbol in coeff_symbols:
                        if coeff_symbol in system_symbols:
                            line += ' ' + coeff_symbol
                
                #print symbols being used with model in appropriate order for symbols
                elif ttype == 'symbols' and (tval is True or tval == 'True'):
                    for system_symbol in system_symbols:
                        if system_symbol in coeff_symbols:
                            line += ' ' + system_symbol
                        else:
                            line += ' NULL'
                            
        return line     