from atomman.models import DataModelDict
import mendeleev
import os

class Potential():
    """class for building LAMMPS input lines from a LAMMPS-potential data model.""" 
    
    def __init__(self, datamodel, pot_dir=None):
        """
        initializes an instance associated with a LAMMPS-potential data model.
        
        Arguments:
        datamodel -- a string or file-like obect of a json/xml data model containing a LAMMPS-potential branch.
        pot_dir -- (optional) the directory location of any artifacts associated with the potential.
        """
        
        self.load(datamodel, pot_dir)
            
    def load(self, datamodel, pot_dir=None):
        """
        loads LAMMPS-potential data model info.
        
        Arguments:
        datamodel -- a string or file-like obect of a json/xml data model containing a LAMMPS-potential branch.
        pot_dir -- (optional) the directory location of any artifacts associated with the potential.
        """

        dm = DataModelDict(datamodel).find('LAMMPS-potential')
        assert len(dm) == 1, 'Exactly one LAMMPS potential implementation must be in data model'
        self.__dm = dm[0]
        
        
        if not isinstance(self.__dm['atom'], list):
            self.__dm['atom'] = [self.__dm['atom']]
        
        for atom in self.__dm['atom']:
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
                atom['mass'] = mendeleev.element(str(atom['element'])).mass

            assert isinstance(atom['mass'], float), 'Mass needs to be a number!'
        
        if pot_dir is not None:
            self.pot_dir = pot_dir
        else:
            self.pot_dir = ''
            
    def __str__(self):
        """The string of the Potential returns its human-readable id"""
        return self.__dm['potential']['id']
    
    @property
    def id(self):
        """The human-readable id for the Potential."""
        return self.__dm['potential']['id']
    
    @property    
    def key(self):
        """The unique hash-key for the Potential."""
        return self.__dm['potential']['key']
    
    @property
    def units(self):
        """The LAMMPS units for the Potential."""
        return self.__dm['units']
    
    @property
    def atom_style(self):
        """The LAMMPS atom_style for the Potential."""
        return self.__dm['atom_style']
    
    @property
    def symbols(self):
        """List of all element-model symbols for the Potential."""
        symbols = []
        atoms = self.__dm['atom']
        for atom in atoms:
            symbols.append(str(atom['symbol']))    
        return symbols
    
    def elements(self, symbols=None):
        """
        Returns the element names associated with a list of element-model symbols.
        
        Arguments:
        symbols -- (optional) list of element-model symbols being used.  If None (default), will use all of the Potential's elements.
        """
        if symbols is None:
            symbols = self.symbols
        
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]
        
        elements = []
        for symbol in symbols:
            for atom in self.__dm['atom']:
                if symbol == atom['symbol']:
                    elements.append(str(atom['element']))
                    break
        assert len(symbols) == len(elements), 'Not all elements found!'
        
        return elements
            
    def masses(self, symbols=None):
        """
        Returns the element masses associated with a list of element-model symbols.
        
        Arguments:
        symbols -- (optional) list of element-model symbols being used.  If None (default), will use all of the Potential's elements.
        """
        if symbols is None:
            symbols = self.symbols
        
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]
        
        masses = []
        for symbol in symbols:
            for atom in self.__dm['atom']:
                if symbol == atom['symbol']:
                    masses.append(atom['mass'])
                    break
        assert len(symbols) == len(masses), 'Not all masses found!'
        
        return masses
           
    def pair_info(self, symbols = None):
        """
        Returns all LAMMPS input command lines associated with the Potential and a list of element-model symbols.
        
        Arguments:
        symbols -- (optional) list of element-model symbols being used.  If None (default), will use all of the Potential's elements.
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
        try:
            style += self.__pair_terms(self.__dm['pair_style']['term']) + '\n'
        except:
            style += '\n'
       
        #Generate pair_coeff lines
        coeff = ''
        coeff_terms = self.__dm['pair_coeff']
        if not isinstance(coeff_terms, list):
            coeff_terms = [coeff_terms]
        for coeff_line in coeff_terms:
            try:
                test = coeff_line['interaction']['symbol']
            except:
                coeff_line['interaction'] = DataModelDict([('symbol', ['*', '*'])])
            
            
            #Always include coeff lines that act on all atoms in the system
            if coeff_line['interaction']['symbol'] == ['*', '*']:
                coeff_symbols = self.symbols
                coeff += 'pair_coeff * *' + self.__pair_terms(coeff_line['term'], symbols, coeff_symbols) + '\n'
                continue
            
            #Many-body potentials will contain a symbols term
            many = False
            for term in coeff_line['term']:
                if 'symbols' in term:
                    many = True
                    break
                                                    
            #Treat as a many-body potential
            if many:
                coeff_symbols = coeff_line['interaction']['symbol']
                coeff += 'pair_coeff * *' + self.__pair_terms(coeff_line['term'], symbols, coeff_symbols) + '\n'
                
            #Treat as pair potential
            else:
                coeff_symbols = coeff_line['interaction']['symbol']
                assert len(coeff_symbols) == 2,     'Pair potential interactions need two listed elements'
                
                #Classic eam style is a special case
                if self.__dm['pair_style']['type'] == 'eam':
                    assert coeff_symbols[0] == coeff_symbols[1], 'Only i==j interactions allowed for eam style'
                    for i in xrange( len(symbols) ):
                        if symbols[i] == coeff_symbols[0]:
                            coeff += 'pair_coeff %i %i' % (i + 1, i + 1) + self.__pair_terms(coeff_line['term'], symbols, coeff_symbols) + '\n'
                
                #All other pair potentials
                else:
                    for i in xrange( len(symbols) ):
                        for j in xrange( i, len(symbols) ):
                            if (symbols[i] == coeff_symbols[0] and symbols[j] == coeff_symbols[1]) or (symbols[i] == coeff_symbols[1] and symbols[j] == coeff_symbols[0]):
                               
                                coeff += 'pair_coeff %i %i' % (i + 1, j + 1) + self.__pair_terms(coeff_line['term'], symbols, coeff_symbols) + '\n'

        #generate additional command lines
        command = ''
        try:
            command_terms = self.__dm['command']
            if not isinstance(command_terms, list):
                command_terms = [command_terms]
        except:
            command_terms = []
        for command_line in command_terms:
            command += self.__pair_terms(command_line['term'], symbols, self.symbols).strip() + '\n'

        return mass + style + coeff + command
    
    def __pair_terms(self, terms, system_symbols = [], coeff_symbols = []):
        """utility function used by self.pair_info() for composing lines from terms"""
        line = ''
        if not isinstance(terms, list): 
            terms = [terms]
        for term in terms:
            for ttype, tval in term.iteritems():
                if ttype == 'option' or ttype == 'parameter': 
                    line += ' ' + str(tval)
                    
                elif ttype == 'file':
                    line += ' ' + str( os.path.join(self.pot_dir, tval) )
                    
                elif ttype == 'symbolsList' and tval == 'True':
                    for coeff_symbol in coeff_symbols:
                        if coeff_symbol in system_symbols:
                            line += ' ' + coeff_symbol
                        
                elif ttype == 'symbols' and tval == 'True':
                    for system_symbol in system_symbols:
                        if system_symbol in coeff_symbols:
                            line += ' ' + system_symbol
                        else:
                            line += ' NULL'
                            
        return line     